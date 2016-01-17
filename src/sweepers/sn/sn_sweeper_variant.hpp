#pragma once

#include "pugixml.hpp"

#include "core/angular_quadrature.hpp"
#include "core/blitz_typedefs.hpp"
#include "core/coarse_data.hpp"
#include "core/core_mesh.hpp"
#include "core/error.hpp"
#include "core/files.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"
#include "core/utils.hpp"
#include "core/xs_mesh_homogenized.hpp"

#include "sweepers/transport_sweeper.hpp"

#include "sn_boundary.hpp"
#include "sn_current_worker.hpp"
#include "sn_sweeper.hpp"

namespace mocc { namespace sn {

    /**
     * The \ref SnSweeperVariant allows for the templating of an Sn sweeper upon
     * a specific differencing scheme. They are derived from the\ref SnSweeper
     * class itself so that the \ref SnSweeperFactory may return a pointer of
     * that type for use elsewhere, so that the template parameter need not be
     * known to the client code. In uses where the differencing scheme is indeed
     * known, the client code may instantiate their own sweeper of this class
     * and have access to a fully-typed \ref CellWorker member. This is useful
     * in the \ref PlaneSweeper_2D3D class, which knows what type of \ref
     * CellWorker it is using.
     */
    template <class Worker>
    class SnSweeperVariant: public SnSweeper {
    public:
        SnSweeperVariant( const pugi::xml_node& input, const CoreMesh& mesh ):
                SnSweeper( input, mesh ),
                cell_worker_( mesh_, ang_quad_ )
        {

            
            return;
        }

        ~SnSweeperVariant() { }

        Worker* worker() {
            return &cell_worker_;
        }

        void sweep( int group ) {
            /// \todo add an is_read() method to the worker, and make sure that
            /// it's ready before continuing.

            // Store the transport cross section somewhere useful
            for( auto &xsr: *xs_mesh_ ) {
                real_t xstr = xsr.xsmactr()[group];
                for( auto &ireg: xsr.reg() ) {
                    xstr_[ireg] = xstr;
                }
            }

            flux_1g_.reference( flux_(blitz::Range::all(), group) );

            // Perform inner iterations
            for( size_t inner=0; inner<n_inner_; inner++ ) {
                // Set the source (add upscatter and divide by 4PI)
                source_->self_scatter( group );
                if( inner == n_inner_-1 && coarse_data_ ) {
                    // Wipe out the existing currents
                    coarse_data_->zero_data( group );
                    this->sweep_1g<sn::Current>( group );
                    coarse_data_->set_has_axial_data(true);
                    coarse_data_->set_has_radial_data(true);
                } else {
                    this->sweep_1g<sn::NoCurrent>( group );
                }
            }

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
        template <typename CurrentWorker>
        void sweep_1g( int group ) {
            CurrentWorker cw( coarse_data_, &mesh_ );
            flux_1g_ = 0.0;
            cell_worker_.set_group( group );

            int nx = mesh_.nx();
            int ny = mesh_.ny();
            int nz = mesh_.nz();

            ArrayF x_flux(ny*nz);
            ArrayF y_flux(nx*nz);
            ArrayF z_flux(nx*ny);

            int iang = 0;
            for( auto ang: ang_quad_ ) {
                // Configure the current worker for this angle
                cw.set_octant( iang / ang_quad_.ndir_oct() + 1 );

                // Get the source for this angle
                auto& q = source_->get_transport( iang );


                cell_worker_.set_angle( iang, ang );

                real_t wgt = ang.weight * HPI;
                real_t ox = ang.ox;
                real_t oy = ang.oy;
                real_t oz = ang.oz;

                // Configure the loop direction. Could template the below for
                // speed at some point.
                int sttx = 0;
                int stpx = nx;
                int xdir = 1;
                if( ox < 0.0 ) {
                    ox = -ox;
                    sttx = nx-1;
                    stpx = -1;
                    xdir = -1;
                }

                int stty = 0;
                int stpy = ny;
                int ydir = 1;
                if( oy < 0.0 ) {
                    oy = -oy;
                    stty = ny-1;
                    stpy = -1;
                    ydir = -1;
                }

                int sttz = 0;
                int stpz = nz;
                int zdir = 1;
                if( oz < 0.0 ) {
                    oz = -oz;
                    sttz = nz-1;
                    stpz = -1;
                    zdir = -1;
                }

                // initialize upwind condition
                x_flux = bc_in_.get_face( group, iang, Normal::X_NORM );
                y_flux = bc_in_.get_face( group, iang, Normal::Y_NORM );
                z_flux = bc_in_.get_face( group, iang, Normal::Z_NORM );

                cw.upwind_work( x_flux, y_flux, z_flux, ang, group);

                for( int iz=sttz; iz!=stpz; iz+=zdir ) {
                    cell_worker_.set_z(iz);
                    for( int iy=stty; iy!=stpy; iy+=ydir ) {
                        cell_worker_.set_y(iy);
                        for( int ix=sttx; ix!=stpx; ix+=xdir ) {
                            // Gross. really need an Sn mesh abstraction
                            real_t psi_x = x_flux[ny*iz + iy];
                            real_t psi_y = y_flux[nx*iz + ix];
                            real_t psi_z = z_flux[nx*iy + ix];

                            size_t i = mesh_.coarse_cell( Position( ix, iy, iz ) );

                            real_t psi = cell_worker_.evaluate( psi_x, psi_y,
                                    psi_z, q[i], xstr_[i], i );
                            //real_t psi = cell_worker_.evaluate_2d( psi_x, psi_y,
                            //        q[i], xstr_[i], i );


                            x_flux[ny*iz + iy] = psi_x;
                            y_flux[nx*iz + ix] = psi_y;
                            z_flux[nx*iy + ix] = psi_z;

                            flux_1g_(i) += psi*wgt;

                            // Stash currents (or not, depending on the
                            // CurrentWorker template parameter)
                            cw.current_work( x_flux[ny*iz + iy],
                                             y_flux[nx*iz + ix],
                                             z_flux[nx*iy + ix],
                                             i, ang, group );
                        }
                    }
                }

                // store the downwind boundary condition
                bc_out_.set_face(0, iang, Normal::X_NORM, x_flux);
                bc_out_.set_face(0, iang, Normal::Y_NORM, y_flux);
                bc_out_.set_face(0, iang, Normal::Z_NORM, z_flux);
                if( gs_boundary_ ) {
                    bc_in_.update( group, iang, bc_out_ );
                }
                iang++;
            }
            // Update the boundary condition
            if( !gs_boundary_ ) {
                bc_in_.update( group, bc_out_ );
            }

            return;
        }


    private:
        Worker cell_worker_;
    };
} }
