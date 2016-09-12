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

#include "util/global_config.hpp"
#include "core/angular_quadrature.hpp"
#include "core/coarse_data.hpp"
#include "core/mesh.hpp"
#include "core/xs_mesh_homogenized.hpp"
#include "sweepers/moc/moc_current_worker.hpp"
#include "sweepers/moc/ray.hpp"
#include "sweepers/moc/ray_data.hpp"
#include "correction_data.hpp"

namespace mocc {
namespace cmdo {
/**
 * See documentation for \ref moc::NoCurrent for canonical documentation
 * for each of the methods.
 */
class CurrentCorrections : public moc::Current {
public:
    CurrentCorrections(CoarseData *coarse_data, const Mesh *mesh,
                       CorrectionData *corrections, const VectorX &qbar,
                       ExpandedXS &xstr_true, ExpandedXS &xstr_split,
                       ExpandedXS &xstr_sn, const AngularQuadrature &ang_quad,
                       const moc::RayData &rays)
        : moc::Current(coarse_data, mesh),
          corrections_(corrections),
          qbar_(qbar),
          xstr_true_(xstr_true),
          xstr_split_(xstr_split),
          xstr_sn_(xstr_sn),
          ang_quad_(ang_quad),
          surf_sum_(mesh_->n_surf_plane() * 2),
          vol_sum_(mesh_->n_cell_plane() * 2),
          vol_norm_(mesh_->n_cell_plane()),
          sigt_sum_(mesh_->n_cell_plane() * 2),
          surf_norm_(mesh_->n_surf_plane() * 2),
          rays_(rays)
    {
        assert(xstr_fm.size() == (int)mesh->n_reg(MeshTreatment::PLANE));
        assert(xstr_sn.size() == (int)mesh->n_reg(MeshTreatment::PIN));
        assert(qbar_.size() == (int)mesh->n_reg(MeshTreatment::PLANE));

        /// \todo this is a mess. come up with a better way to interact with the
        /// cross sections or make this structure available through the mesh or
        /// something
        mplane_offset_.reserve(mesh_->macroplane_index().back());
        int mplane          = 0;
        int offset          = 0;
        int cells_per_plane = mesh_->nx() * mesh_->ny();
        mplane_offset_.push_back(0);
        for (const auto index : mesh->macroplane_index()) {
            if (index != mplane) {
                mplane_offset_.push_back(offset);
                mplane = index;
            }
            offset += cells_per_plane;
        }
        return;
    }

    inline void set_group(int group)
    {
        group_    = group;
        residual_ = {{0.0, 0.0, 0.0}};
    }

    std::array<real_t, 3> residual() const
    {
        auto resid = residual_;
        resid[0]   = std::sqrt(resid[0]);
        resid[1]   = std::sqrt(resid[1]);
        resid[2]   = std::sqrt(resid[2]);
        return resid;
    };

    /**
     * \brief Set up the worker to treat the indexed macroplane
     *
     * \param plane the index of the macroplane that is about to be swept
     */
    MOCC_FORCE_INLINE void set_plane(int plane)
    {
        assert(plane < (int)mplane_offset_.size());
        Current::set_plane(plane);
        cell_offset_xs_ = mplane_offset_[plane];

        return;
    }

    inline void post_ray(const FluxStore &psi1, const FluxStore &psi2,
                         const ArrayB1 &e_tau, const moc::Ray &ray,
                         int first_reg)
    {
#pragma omp critical
        {
            auto all          = blitz::Range::all();
            auto current      = coarse_data_->current(all, group_);
            auto surface_flux = coarse_data_->surface_flux(all, group_);

            int cell_fw = ray.cm_cell_fw();
            int cell_bw = ray.cm_cell_bw();
            int surf_fw = ray.cm_surf_fw();
            int surf_bw = ray.cm_surf_bw();
            int iseg_fw = 0;
            int iseg_bw = ray.nseg();

            // Use an offset for the current contributions, but NOT for the
            // correction factors, which use plane-by-plane indexing
            int norm_fw = (int)mesh_->surface_normal(surf_fw);
            int norm_bw = (int)mesh_->surface_normal(surf_bw);
            current(surf_fw + surf_offset_) +=
                psi1[iseg_fw] * current_weights_[norm_fw];
            current(surf_bw + surf_offset_) -=
                psi2[iseg_bw] * current_weights_[norm_bw];
            surface_flux(surf_fw + surf_offset_) +=
                psi1[iseg_fw] * flux_weights_[norm_fw];
            surface_flux(surf_bw + surf_offset_) -=
                psi2[iseg_bw] * flux_weights_[norm_bw];

            surf_sum_(surf_fw * 2 + 0) += psi1[iseg_fw];
            surf_sum_(surf_bw * 2 + 1) += psi2[iseg_bw];
            surf_norm_(surf_fw * 2 + 0) += 1.0;
            surf_norm_(surf_bw * 2 + 1) += 1.0;

            auto begin = ray.cm_data().cbegin();
            auto end   = ray.cm_data().cend();
            for (auto crd = begin; crd != end; ++crd) {
                // Hopefully branch prediction saves me here.
                if (crd->fw != Surface::INVALID) {
                    // Store forward volumetric stuff
                    for (unsigned i = 0; i < crd->nseg_fw; i++) {
                        int ireg         = ray.seg_index(iseg_fw) + first_reg;
                        real_t xstr      = xstr_split_[ireg];
                        real_t xstr_true = xstr_true_[ireg];
                        real_t t = ang_.rsintheta * ray.seg_len(iseg_fw);
                        real_t fluxvol =
                            t * qbar_(ireg) +
                            (psi1[iseg_fw] - psi1[iseg_fw + 1]) / xstr;
                        if (fluxvol < 0.0) {
                            std::cout << "negative psi-bar: " << ireg << " "
                                      << iseg_fw << " " << fluxvol << "\n";
                            std::cout << t << " " << qbar_(ireg) << " "
                                      << psi1[iseg_fw] << " "
                                      << psi1[iseg_fw + 1] << " " << xstr << " "
                                      << e_tau(iseg_fw) << " "
                                      << 1.0 - std::exp(-xstr * t) << "\n";
                            throw EXCEPT("neg psibar");
                        }
                        vol_sum_(cell_fw * 2 + 0) += fluxvol;
                        vol_norm_(cell_fw) += t;
                        sigt_sum_(cell_fw * 2 + 0) += xstr_true * fluxvol;
                        iseg_fw++;
                    }
                    // Store FW surface stuff
                    norm_fw = (int)surface_to_normal(crd->fw);
                    surf_fw = mesh_->coarse_surf(cell_fw, crd->fw);
                    current(surf_fw + surf_offset_) +=
                        psi1[iseg_fw] * current_weights_[norm_fw];
                    surface_flux(surf_fw + surf_offset_) +=
                        psi1[iseg_fw] * flux_weights_[norm_fw];
                    surf_sum_(surf_fw * 2 + 0) += psi1[iseg_fw];
                    surf_norm_(surf_fw * 2 + 0) += 1.0;
                }

                if (crd->bw != Surface::INVALID) {
                    // Store backward volumetric stuff
                    for (unsigned i = 0; i < crd->nseg_bw; i++) {
                        iseg_bw--;
                        int ireg         = ray.seg_index(iseg_bw) + first_reg;
                        real_t xstr      = xstr_split_[ireg];
                        real_t xstr_true = xstr_true_[ireg];
                        real_t t       = ang_.rsintheta * ray.seg_len(iseg_bw);
                        real_t fluxvol = t * qbar_(ireg) +
                                         e_tau(iseg_bw) *
                                             (psi2[iseg_bw + 1] - qbar_(ireg)) /
                                             xstr;
                        vol_sum_(cell_bw * 2 + 1) += fluxvol;
                        sigt_sum_(cell_bw * 2 + 1) += xstr_true * fluxvol;
                    }
                    // Store BW surface stuff
                    norm_bw = (int)surface_to_normal(crd->bw);
                    surf_bw = mesh_->coarse_surf(cell_bw, crd->bw);
                    current(surf_bw + surf_offset_) -=
                        psi2[iseg_bw] * current_weights_[norm_bw];
                    surface_flux(surf_bw + surf_offset_) -=
                        psi2[iseg_bw] * flux_weights_[norm_bw];
                    surf_sum_(surf_bw * 2 + 1) += psi2[iseg_bw];
                    surf_norm_(surf_bw * 2 + 1) += 1.0;
                }

                cell_fw = mesh_->coarse_neighbor(cell_fw, (crd)->fw);
                cell_bw = mesh_->coarse_neighbor(cell_bw, (crd)->bw);
            }
        }
        return;
    }

    inline void set_angle(Angle ang, real_t spacing)
    {
        moc::Current::set_angle(ang, spacing);
        ang_ = ang;

        // Zero out all of the flux sum arrays
        surf_sum_  = 0.0;
        surf_norm_ = 0.0;
        vol_sum_   = 0.0;
        vol_norm_  = 0.0;
        sigt_sum_  = 0.0;

#pragma omp barrier
        return;
    }

    void post_angle(int iang)
    {
#pragma omp single
        {
            // Do the stock area normailzation
            moc::Current::post_angle(iang);

            // Normalize the flux and sigt values and calculate
            // correction factors for the current angle/energy
            for (size_t i = 0; i < vol_norm_.size(); i++) {
                sigt_sum_(2 * i + 0) /= vol_sum_(2 * i + 0);
                sigt_sum_(2 * i + 1) /= vol_sum_(2 * i + 1);
                vol_sum_(2 * i + 0) /= vol_norm_(i);
                vol_sum_(2 * i + 1) /= vol_norm_(i);
            }

            // Doing area-based normalization
            // surf_sum_ /= surf_norm_;

            calculate_corrections(iang, group_);
        }
#pragma omp barrier
        return;
    }

private:
    CorrectionData *corrections_;
    // References to the source and cross sections as defined on the fine mesh.
    // We need these to get actual angular flux for a ray segment
    const VectorX &qbar_;
    // The actual fine-mesh cross sections. These should not include any
    // modifications due to TL splitting or similar. Specifically, they should
    // be consistent with what would be used  to perform scalar flux-weighed
    // cross sections.
    ExpandedXS xstr_true_;
    // The fine-mesh cross sections that are used in the associated MoC sweep
    // procedure. These should include all modifications that the MoC sweeper
    // uses. Specifically, these should be consistent with the cross sections
    // used in the MoC sweep kernel to propagate flux, since we use them to
    // convert from Delta Psi to Psi bar
    ExpandedXS xstr_split_;

    ExpandedXS xstr_sn_;

    unsigned cell_offset_xs_;

    const AngularQuadrature &ang_quad_;

    ArrayB1 surf_sum_;
    ArrayB1 vol_sum_;
    ArrayB1 vol_norm_;
    ArrayB1 sigt_sum_;
    ArrayB1 surf_norm_;

    Angle ang_;

    const moc::RayData &rays_;

    std::array<real_t, 3> residual_;

    // Index offset to get the first coarse mesh region in a given macroplane
    std::vector<int> mplane_offset_;

    /** \page surface_norm Surface Normalization
     * Surface normalization \todo discuss surface normalization
     */
    void calculate_corrections(size_t ang, size_t group);
};
}
}
