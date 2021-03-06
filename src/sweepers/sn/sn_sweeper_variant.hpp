/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include "pugixml.hpp"
#include "util/blitz_typedefs.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/global_config.hpp"
#include "util/utils.hpp"
#include "core/angular_quadrature.hpp"
#include "core/coarse_data.hpp"
#include "core/core_mesh.hpp"
#include "core/mesh.hpp"
#include "core/transport_sweeper.hpp"
#include "core/xs_mesh_homogenized.hpp"
#include "sn_current_worker.hpp"
#include "sn_sweeper.hpp"

namespace mocc {
namespace sn {
/**
 * The \ref SnSweeperVariant allows for the templating of an Sn sweeper upon
 * a specific differencing scheme. They are derived from the \ref SnSweeper
 * class itself so that the \ref SnSweeperFactory may return a pointer of
 * that type for use elsewhere, so that the template parameter need not be
 * known to the client code. In uses where the differencing scheme is indeed
 * known, the client code may instantiate their own sweeper of this class
 * and have access to a fully-typed \ref CellWorker member. This is useful
 * in the \ref cmdo::PlaneSweeper_2D3D class, which knows what type of \ref
 * CellWorker it is using.
 *
 * Specialization is carried out using the Curiously-Recurring Template pattern.
 * See the specific implementations (e.g. \ref SnSweeper_DD or \ref
 * cmdo::SnSweeper_CDD).
 */
template <class Equation> class SnSweeperVariant : public SnSweeper {
public:
    /**
     * \todo Use some more template magic to make the members of the TreadState
     * sweeper-dependent
     */
    struct ThreadState {
    public:
        real_t ty;
        real_t tz;
        int iang;
        int iang_2d;
        int macroplane;
        real_t ox;
        real_t oy;
        real_t oz;
        Angle angle;
    };
    SnSweeperVariant(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeper(input, mesh), plane_size_(mesh.nx() * mesh_.ny())
    {
        return;
    }

    ~SnSweeperVariant()
    {
    }

    inline real_t evaluate_2d(real_t &flux_x, real_t &flux_y, real_t q,
                              real_t xstr, int i,
                              const ThreadState &t_state) const
    {
        return static_cast<const Equation &>(*this).evaluate_2d(
            flux_x, flux_y, q, xstr, i, t_state);
    }
    inline real_t evaluate(real_t &flux_x, real_t &flux_y, real_t &flux_z,
                           real_t q, real_t xstr, int i,
                           const ThreadState &t_state) const
    {
        return static_cast<const Equation &>(*this).evaluate(
            flux_x, flux_y, flux_z, q, xstr, i, t_state);
    }

    void sweep(int group) override
    {
        assert(source_);
        timer_.tic();

        timer_xsupdate_.tic();
        xs_mesh_->update();
        timer_xsupdate_.toc();

        timer_sweep_.tic();

        group_ = group;

        /// \todo add an is_ready() method to the worker, and make sure that
        /// it's ready before continuing.

        // Store the transport cross section somewhere useful
        xstr_.expand(group);

        flux_1g_.reference(flux_(blitz::Range::all(), group));

        // Perform inner iterations
        for (unsigned inner = 0; inner < n_inner_; inner++) {
            // Set the source (add upscatter and divide by 4PI)
            source_->self_scatter(group);
            if (inner == n_inner_ - 1 && coarse_data_) {
                // Wipe out the existing currents
                coarse_data_->zero_data(group);
                coarse_data_->source() = "Sn Sweeper";
                if (core_mesh_->is_2d()) {
                    this->sweep_1g_2d<sn::Current>(group);
                } else {
                    this->sweep_1g<sn::Current>(group);
                    coarse_data_->set_has_axial_data(true);
                }
                coarse_data_->set_has_radial_data(true);
            } else {
                if (core_mesh_->is_2d()) {
                    this->sweep_1g_2d<sn::NoCurrent>(group);
                } else {
                    this->sweep_1g<sn::NoCurrent>(group);
                }
            }
        }

        // Clean up zeros
            int n_neg = 0;
        for(auto &v: flux_1g_){
            if(v < 0) {
                ++n_neg;
                v = 0.0;
            }
        }
        if(n_neg > 0) {
            LogFile << "Fixed " << n_neg << "negative Sn fluxes\n";
        }

        timer_.toc();
        timer_sweep_.toc();
        return;
    }

protected:
    /**
     * \brief Generic Sn sweep procedure for orthogonal mesh.
     *
     * This routine performs a single, one-group transport sweep with Sn. It
     * is templated on two parameters to tailor it to different differencing
     * schemes and current calculation requirements. To see examples of
     * template parameters, look at \ref sn::CellWorker and \ref
     * sn::Current.
     */
    template <typename CurrentWorker> void sweep_1g(int group)
    {
        flux_1g_ = 0.0;
#pragma omp parallel default(shared)
        {
            CurrentWorker cw(coarse_data_, &mesh_);

            ArrayB1 t_flux(n_reg_);
            t_flux = 0.0;

            int nx = mesh_.nx();
            int ny = mesh_.ny();
            int nz = mesh_.nz();

            real_t *x_flux;
            real_t *y_flux;
            real_t *z_flux;

            Angle angle;
            ThreadState t_state;

#pragma omp for
            for (int iang = 0; iang < ang_quad_.ndir(); iang++) {
                angle           = ang_quad_[iang];
                t_state.iang    = iang;
                t_state.iang_2d = iang % (ang_quad_.ndir() / 2);
                t_state.angle   = ang_quad_[iang];
                // Configure the current worker for this angle
                cw.set_octant(angle);

                // Get the source for this angle
                auto &q = source_->get_transport(iang);

                real_t wgt = angle.weight * HPI;
                t_state.ox = angle.ox;
                t_state.oy = angle.oy;
                t_state.oz = angle.oz;

                // Configure the loop direction. Could template the below for
                // speed at some point.
                int sttx = 0;
                int stpx = nx;
                int xdir = 1;
                if (t_state.ox < 0.0) {
                    t_state.ox = -t_state.ox;
                    sttx       = nx - 1;
                    stpx       = -1;
                    xdir       = -1;
                }

                int stty = 0;
                int stpy = ny;
                int ydir = 1;
                if (t_state.oy < 0.0) {
                    t_state.oy = -t_state.oy;
                    stty       = ny - 1;
                    stpy       = -1;
                    ydir       = -1;
                }

                int sttz = 0;
                int stpz = nz;
                int zdir = 1;
                if (t_state.oz < 0.0) {
                    t_state.oz = -t_state.oz;
                    sttz       = nz - 1;
                    stpz       = -1;
                    zdir       = -1;
                }

                // initialize upwind condition
                auto xf = bc_out_.get_face(0, iang, Normal::X_NORM);
                x_flux  = xf.second;
                auto yf = bc_out_.get_face(0, iang, Normal::Y_NORM);
                y_flux  = yf.second;
                auto zf = bc_out_.get_face(0, iang, Normal::Z_NORM);
                z_flux  = zf.second;
                bc_in_.copy_face(group, iang, Normal::X_NORM, x_flux);
                bc_in_.copy_face(group, iang, Normal::Y_NORM, y_flux);
                bc_in_.copy_face(group, iang, Normal::Z_NORM, z_flux);

                cw.upwind_work(x_flux, y_flux, z_flux, angle, group);

                for (int iz = sttz; iz != stpz; iz += zdir) {
                    t_state.tz         = t_state.oz / mesh_.dz(iz);
                    t_state.macroplane = macroplanes_[iz];
                    for (int iy = stty; iy != stpy; iy += ydir) {
                        t_state.ty = t_state.oy / mesh_.dy(iy);
                        for (int ix = sttx; ix != stpx; ix += xdir) {
                            // Gross. really need an Sn mesh abstraction
                            real_t psi_x = x_flux[ny * iz + iy];
                            real_t psi_y = y_flux[nx * iz + ix];
                            real_t psi_z = z_flux[nx * iy + ix];

                            int i = mesh_.coarse_cell(Position(ix, iy, iz));

                            real_t psi =
                                this->evaluate(psi_x, psi_y, psi_z, q[i],
                                               xstr_[i], i, t_state);

                            x_flux[ny * iz + iy] = psi_x;
                            y_flux[nx * iz + ix] = psi_y;
                            z_flux[nx * iy + ix] = psi_z;

                            t_flux(i) += psi * wgt;

                            cw.current_work(
                                x_flux[ny * iz + iy], y_flux[nx * iz + ix],
                                z_flux[nx * iy + ix], i, angle, group);
                        }
                    }
                }

                if (gs_boundary_) {
                    bc_in_.update(group, iang, bc_out_);
                }
            } // Angles
              // Update the boundary condition
#pragma omp single
            if (!gs_boundary_) {
                bc_in_.update(group, bc_out_);
            }

// Reduce scalar flux
#pragma omp critical
            {
                flux_1g_ += t_flux;
            }
        } // OMP Parallel

        return;
    } // sweep_1g (3-D)

    /**
     * \brief Generic Sn sweep procedure for 2-D orthogonal mesh.
     *
     * This routine performs a single, one-group transport sweep with Sn. It
     * is templated on two parameters to tailor it to different differencing
     * schemes and current calculation requirements. To see examples of
     * template parameters, look at \ref sn::CellWorker and \ref
     * sn::Current.
     */
    template <typename CurrentWorker> void sweep_1g_2d(int group)
    {
        flux_1g_ = 0.0;
#pragma omp parallel default(shared)
        {
            CurrentWorker cw(coarse_data_, &mesh_);

            ArrayB1 t_flux(n_reg_);
            t_flux = 0.0;

            int nx = mesh_.nx();
            int ny = mesh_.ny();

            real_t *x_flux;
            real_t *y_flux;

            Angle angle;
            ThreadState t_state;

            t_state.macroplane = 0;

#pragma omp for
            for (int iang = 0; iang < ang_quad_.ndir() / 2; iang++) {
                angle           = ang_quad_[iang];
                t_state.iang    = iang;
                t_state.angle   = angle;
                t_state.iang_2d = iang % (ang_quad_.ndir() / 2);
                // Configure the current worker for this angle
                cw.set_octant(angle);

                // Get the source for this angle
                auto &q = source_->get_transport(iang);

                real_t wgt = angle.weight * PI;
                t_state.ox = angle.ox;
                t_state.oy = angle.oy;
                t_state.oz = angle.oz;

                // Configure the loop direction. Could template the below for
                // speed at some point.
                int sttx = 0;
                int stpx = nx;
                int xdir = 1;
                if (t_state.ox < 0.0) {
                    t_state.ox = -t_state.ox;
                    sttx       = nx - 1;
                    stpx       = -1;
                    xdir       = -1;
                }

                int stty = 0;
                int stpy = ny;
                int ydir = 1;
                if (t_state.oy < 0.0) {
                    t_state.oy = -t_state.oy;
                    stty       = ny - 1;
                    stpy       = -1;
                    ydir       = -1;
                }

                // initialize upwind condition
                auto xf = bc_out_.get_face(0, iang, Normal::X_NORM);
                x_flux  = xf.second;
                auto yf = bc_out_.get_face(0, iang, Normal::Y_NORM);
                y_flux  = yf.second;
                bc_in_.copy_face(group, iang, Normal::X_NORM, x_flux);
                bc_in_.copy_face(group, iang, Normal::Y_NORM, y_flux);

                cw.upwind_work(x_flux, y_flux, angle, group);

                t_state.tz = t_state.oz / mesh_.dz(0);
                for (int iy = stty; iy != stpy; iy += ydir) {
                    t_state.ty = t_state.oy / mesh_.dy(iy);
                    for (int ix = sttx; ix != stpx; ix += xdir) {
                        // Gross. really need an Sn mesh abstraction
                        real_t psi_x = x_flux[iy];
                        real_t psi_y = y_flux[ix];

                        int i = mesh_.coarse_cell(Position(ix, iy, 0));

                        real_t psi = this->evaluate_2d(psi_x, psi_y, q[i],
                                                       xstr_[i], i, t_state);

                        x_flux[iy] = psi_x;
                        y_flux[ix] = psi_y;

                        t_flux(i) += psi * wgt;

                        cw.current_work(x_flux[iy], y_flux[ix], i, angle,
                                        group);
                    }
                }

                if (gs_boundary_) {
                    bc_in_.update(group, iang, bc_out_);
                }
            } // Angles
              // Update the boundary condition
#pragma omp single
            if (!gs_boundary_) {
                bc_in_.update(group, bc_out_);
            }

// Reduce scalar flux
#pragma omp critical
            {
                flux_1g_ += t_flux;
            }
        } // OMP Parallel

        return;
    }

    int plane_size_;
    int group_;
};
}
}
