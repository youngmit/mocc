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

#include <cmath>
#include "util/force_inline.hpp"
#include "util/global_config.hpp"
#include "core/constants.hpp"
#include "core/geometry/angle.hpp"
#include "ray.hpp"

namespace mocc {
namespace moc {
/**
 * This can be used as a template parameter to the \ref
 * MoCSweeper::sweep1g() method. Using this class in such a way avoids
 * the extra work needed to compute currents, and with any optimization
 * enabled, should yield code identical to a hand-written MoC sweep
 * without the current work.
 */
class NoCurrent {
public:
    /**
     * \brief Subscriptable abstraction for only storing a scalar value
     *
     * This class allows the MoC sweeper kernel to be agnostic to the type of
     * storage needed to represent the flux along a ray. In cases where current
     * or some other value is needed from the sweeper, it is necessary to keep
     * the angular flux along the entire length of the ray. In other situations
     * where this is unnecessary, it is a waste to keep track of this ray flux,
     * and sufficient to just maintain the angular flux at the furthest-swept
     * position on the ray. To allow the sweeper kernel to be written in a
     * manner allowing both options, this class implements a subscript operator,
     * which points to the same scalar every time, which should be elided by an
     * optimizing compiler.
     *
     * \sa moc::Current::FluxStore
     */
    class FluxStore {
    public:
        FluxStore(int size)
        {
            return;
        }
        real_t &operator[](int i)
        {
            return psi_;
        }
        real_t operator[](int i) const
        {
            return psi_;
        }

    private:
        real_t psi_;
    };

    NoCurrent()
    {
    }
    NoCurrent(CoarseData *data, const Mesh *mesh)
    {
    }

    /**
     * Defines work to be done following the sweep of a single ray. This
     * is useful for when you need to do something with the angular
     * flux.
     */
    MOCC_FORCE_INLINE void post_ray(FluxStore psi1, FluxStore psi2,
                                    const ArrayB1 &e_tau, const Ray &ray,
                                    int first_reg)
    {
        return;
    }

    /**
     * Defines work to be done before sweeping rays in a given angle.
     */
    MOCC_FORCE_INLINE void set_angle(Angle ang, real_t spacing)
    {
        return;
    }

    /**
     * Defines work to be done after sweeping all rays in a given angle.
     */
    MOCC_FORCE_INLINE void post_angle(int iang)
    {
        return;
    }

    MOCC_FORCE_INLINE void set_plane(int iplane)
    {
        return;
    }

    MOCC_FORCE_INLINE void post_sweep()
    {
        return;
    }

    MOCC_FORCE_INLINE void post_plane()
    {
        return;
    }

    MOCC_FORCE_INLINE void set_group(int group)
    {
        return;
    }
};

/**
 * This class can be used as a template parameter to the \ref
 * MoCSweeper::sweep1g() method to control whether or not extra work is
 * done during the sweep to compute currents. Specifically, when this
 * class is used as the template parameter, currents are calculated.
 *
 * See documentation for \ref moc::NoCurrent for canonical documentation
 * for each of the methods.
 */
class Current {
public:
    /**
     * \brief Typedef for a STL vector of real_t for storing angular flux along
     * a ray.
     *
     * This type is used to store the flux along the entire length of the ray
     * when such information is needed from the MoC sweeper kernel.
     *
     * \sa moc::NoCurrent::FluxStore
     */
    typedef std::vector<real_t> FluxStore;

    Current()
        : coarse_data_(nullptr), mesh_(nullptr)
    {
    }

    Current(CoarseData *data, const Mesh *mesh)
        : coarse_data_(data), mesh_(mesh)
    {
        return;
    }

    MOCC_FORCE_INLINE void post_angle(int iang)
    {
        return;
    };

    MOCC_FORCE_INLINE void post_plane()
    {
        return;
    }

    MOCC_FORCE_INLINE void set_group(int group)
    {
        group_ = group;
    }

    MOCC_FORCE_INLINE void set_plane(int plane)
    {
        plane_       = plane;
        cell_offset_ = mesh_->coarse_cell_offset(plane);
        surf_offset_ = mesh_->coarse_surf_offset(plane);
    }

    MOCC_FORCE_INLINE void set_angle(Angle ang, real_t spacing)
    {
#pragma omp single
        {
            // Scale the angle weight to sum to 4*PI
            real_t w = ang.weight * PI;
            // Multiply by dz so that we conform to the actual coarse
            // mesh area.
            real_t dz = mesh_->dz(plane_);

            current_weights_[0] =
                w * ang.ox * spacing / std::abs(std::cos(ang.alpha)) * dz;
            current_weights_[1] =
                w * ang.oy * spacing / std::abs(std::sin(ang.alpha)) * dz;
            flux_weights_[0] = w * spacing / std::abs(std::cos(ang.alpha)) * dz;
            flux_weights_[1] = w * spacing / std::abs(std::sin(ang.alpha)) * dz;
        }
#pragma omp barrier
    }

    MOCC_FORCE_INLINE void post_ray(const FluxStore &psi1,
                                    const FluxStore &psi2, const ArrayB1 &e_tau,
                                    const Ray &ray, int first_reg)
    {
#pragma omp critical
        {
            /**
             * \todo this is going to perform poorly as implemented. this is a
             * really large critical section, which we might be able to do as
             * atomic updates, as is done in the Sn current worker. The Blitz
             * array slicing is not thread safe, though, so we would need to be
             * careful. It'd be nice to fiddle with this and profile once things
             * settle down some.
             */
            auto all          = blitz::Range::all();
            auto current      = coarse_data_->current(all, group_);
            auto surface_flux = coarse_data_->surface_flux(all, group_);

            size_t cell_fw = ray.cm_cell_fw() + cell_offset_;
            size_t cell_bw = ray.cm_cell_bw() + cell_offset_;

            int surf_fw = ray.cm_surf_fw() + surf_offset_;
            int surf_bw = ray.cm_surf_bw() + surf_offset_;
            int iseg_fw = 0;
            int iseg_bw = ray.nseg();

            int norm_fw = (int)mesh_->surface_normal(surf_fw);
            int norm_bw = (int)mesh_->surface_normal(surf_bw);
            current(surf_fw) += psi1[iseg_fw] * current_weights_[norm_fw];
            current(surf_bw) -= psi2[iseg_bw] * current_weights_[norm_bw];
            surface_flux(surf_fw) += psi1[iseg_fw] * flux_weights_[norm_fw];
            surface_flux(surf_bw) += psi2[iseg_bw] * flux_weights_[norm_bw];

            auto begin = ray.cm_data().cbegin();
            auto end   = ray.cm_data().cend();
            for (auto crd = begin; crd != end; ++crd) {
                // Hopefully branch prediction saves me here.
                if (crd->fw != Surface::INVALID) {
                    iseg_fw += crd->nseg_fw;
                    norm_fw = (int)surface_to_normal(crd->fw);
                    surf_fw = mesh_->coarse_surf(cell_fw, crd->fw);
                    current(surf_fw) +=
                        psi1[iseg_fw] * current_weights_[norm_fw];
                    surface_flux(surf_fw) +=
                        psi1[iseg_fw] * flux_weights_[norm_fw];
                }

                if (crd->bw != Surface::INVALID) {
                    iseg_bw -= crd->nseg_bw;
                    norm_bw = (int)surface_to_normal(crd->bw);
                    surf_bw = mesh_->coarse_surf(cell_bw, crd->bw);
                    current(surf_bw) -=
                        psi2[iseg_bw] * current_weights_[norm_bw];
                    surface_flux(surf_bw) +=
                        psi2[iseg_bw] * flux_weights_[norm_bw];
                }

                cell_fw = mesh_->coarse_neighbor(cell_fw, (crd)->fw);
                cell_bw = mesh_->coarse_neighbor(cell_bw, (crd)->bw);
            }
        } // OMP critical
        return;
    }

    /**
     * \brief Clean up anything that needs to be done after sweeping all angles
     *
     * In the context of the \ref CurrentWorker and most of its children, this
     * only includes expanding the currents to the full PIN grid from the
     * potentially smaller MoC axial grid.
     */
    MOCC_FORCE_INLINE void post_sweep()
    {
#pragma omp single
        {
            // Check to see if we need to expand the currents across the mesh.
            if ((int)mesh_->nz() - 1 != (mesh_->macroplane_index().back())) {
                // In the presence of subplaning, the currents coming from the
                // sweeper are stored by macroplane, packed towards the bottom
                // of the mesh.  To safely perform an in-place expansion, we
                // will expand the currents in reverse, filling from the top
                // down. This prevents over-writing of the source currents from
                // the MoC sweep before having a chance to expand them, as would
                // happen if the expansion went from the bottom up.
                int iz = mesh_->nz() - 1;
                for (auto mplane_it = mesh_->macroplane_index().crbegin();
                     mplane_it != mesh_->macroplane_index().crend();
                     ++mplane_it) {
                    int stt_out = mesh_->plane_surf_xy_begin(iz);
                    int stp_out = mesh_->plane_surf_end(iz);
                    int ip      = *mplane_it;
                    int stt_in  = mesh_->plane_surf_xy_begin(ip);
                    int stp_in  = mesh_->plane_surf_end(ip);

                    coarse_data_->current(blitz::Range(stt_out, stp_out),
                                          group_) =
                        coarse_data_->current(blitz::Range(stt_in, stp_in),
                                              group_);

                    iz--;
                }
            }

            auto all          = blitz::Range::all();
            auto current      = coarse_data_->current(all, group_);
            auto surface_flux = coarse_data_->surface_flux(all, group_);
            // Normalize the surface currents
            for (size_t plane = 0; plane < mesh_->nz(); plane++) {
                for (int surf = mesh_->plane_surf_xy_begin(plane);
                     surf != (int)mesh_->plane_surf_end(plane); ++surf) {
                    real_t area = mesh_->coarse_area(surf);
                    current(surf) /= area;
                    surface_flux(surf) /= area;
                }
            }
        }
        return;
    }

protected:
    CoarseData *coarse_data_;
    const Mesh *mesh_;
    std::array<real_t, 2> current_weights_;
    std::array<real_t, 2> flux_weights_;

    int plane_;
    int group_;
    int cell_offset_;
    int surf_offset_;
};
}
}
