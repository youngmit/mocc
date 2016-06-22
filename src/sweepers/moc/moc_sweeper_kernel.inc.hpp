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

// These inclusions are unnecesary, and only provided so Doxygen can chase
// references.


/**
 * \file
 * This contains the actual MoC sweeper kernel, and --fortune willing--
 * the only one.
 */

/**
 * \brief Perform an MoC sweep
 *
 * This method performs a single source iteration using MoC, sweeping
 * all angles and rays once for a given angle. The \c cw parameter
 * allows for auxiliary work to be done during the sweep, without
 * affecting runtime performance when not needed. Examples of this are
 * currents for CMFD coupling and correction factors for 2D3D/CDD
 * coupling. See \ref moc::Current and \ref cmdo::CurrentCorrections
 * for examples of these.
 */
template <typename CurrentWorker> void sweep1g(int group, CurrentWorker &cw)
{
    flux_1g_ = 0.0;

    cw.set_group(group);

#pragma omp parallel default(shared)
    {
        ArrayB1 e_tau(rays_.max_segments());
        ArrayB1 psi1(rays_.max_segments() + 1);
        ArrayB1 psi2(rays_.max_segments() + 1);
        ArrayB1 t_flux(n_reg_);
        t_flux = 0.0;

        int iplane = 0;
        for (const auto plane_ray_id : mesh_.unique_planes()) {
            auto &boundary_in  = boundary_[iplane];
            auto &boundary_out = boundary_out_[iplane];
            cw.set_plane(iplane);
            const auto &plane_rays = rays_[plane_ray_id];
            int first_reg          = mesh_.first_reg_plane(iplane);
            int iang               = 0;
            // Angles
            for (const auto &ang_rays : plane_rays) {
                // Get the source for this angle
                auto &qbar = source_->get_transport(iang);
                int iang1  = iang;
                int iang2  = ang_quad_.reverse(iang);
                Angle ang  = ang_quad_[iang];

                // Get the boundary condition storage
                const real_t *bc_in_1 =
                    boundary_in.get_boundary(group, iang1).second;
                real_t *bc_out_1 = boundary_out.get_boundary(0, iang1).second;
                const real_t *bc_in_2 =
                    boundary_in.get_boundary(group, iang2).second;
                real_t *bc_out_2 = boundary_out.get_boundary(0, iang2).second;

                // Set up the current worker for sweeping this angle
                cw.set_angle(ang, rays_.spacing(iang));

                real_t stheta  = sin(ang.theta);
                real_t rstheta = 1.0 / stheta;
                real_t wt_v_st = ang.weight * rays_.spacing(iang) * stheta * PI;

#pragma omp for schedule(static, 1)
                for (int iray = 0; iray < (int)ang_rays.size(); iray++) {
                    const auto &ray = ang_rays[iray];

                    int bc1 = ray.bc(0);
                    int bc2 = ray.bc(1);
                    assert(bc1 < boundary_in.get_boundary(group, iang1).first);
                    assert(bc2 < boundary_in.get_boundary(group, iang1).first);

                    // Compute exponentials
                    for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                        int ireg    = ray.seg_index(iseg) + first_reg;
                        e_tau(iseg) = 1.0 -
                                      exp_.exp(-xstr_(ireg) *
                                               ray.seg_len(iseg) * rstheta);
                    }

                    // Forward direction
                    // Initialize from bc
                    psi1(0) = bc_in_1[bc1];

                    // Propagate through core geometry
                    for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                        int ireg = ray.seg_index(iseg) + first_reg;
                        real_t psi_diff =
                            (psi1(iseg) - qbar[ireg]) * e_tau(iseg);
                        psi1(iseg + 1) = psi1(iseg) - psi_diff;
                        t_flux(ireg) += psi_diff * wt_v_st;
                    }
                    // Store boundary condition
                    bc_out_1[bc2] = psi1(ray.nseg());

                    // Backward direction
                    // Initialize from bc
                    psi2(ray.nseg()) = bc_in_2[bc2];

                    // Propagate through core geometry
                    for (int iseg = ray.nseg() - 1; iseg >= 0; iseg--) {
                        int ireg = ray.seg_index(iseg) + first_reg;
                        real_t psi_diff =
                            (psi2(iseg + 1) - qbar[ireg]) * e_tau(iseg);
                        psi2(iseg) = psi2(iseg + 1) - psi_diff;
                        t_flux(ireg) += psi_diff * wt_v_st;
                    }
                    // Store boundary condition
                    bc_out_2[bc1] = psi2(0);

                    // Stash currents
                    cw.post_ray(psi1, psi2, e_tau, ray, first_reg);
                } // Rays
                cw.post_angle(iang);

                // Try tasks?
                if (gauss_seidel_boundary_)
#pragma omp single
                {
                    boundary_in.update(group, iang1, boundary_out);
                    boundary_in.update(group, iang2, boundary_out);
                }

                iang++;
            } // angles
            if (!gauss_seidel_boundary_)
#pragma omp single
            {
                boundary_in.update(group, boundary_out);
            }

            iplane++;
        } // planes

#pragma omp barrier
#pragma omp critical
        {
            /// \todo make an array operation after refactoring out valarray
            for (int i = 0; i < (int)n_reg_; i++) {
                flux_1g_(i) += t_flux(i);
            }
        }
#pragma omp barrier
// Scale the scalar flux by the volume and add back the source
#pragma omp single
        {
            // \todo this is not correct for angle-dependent sources!
            auto &qbar = source_->get_transport(0);
            for (int i = 0; i < (int)n_reg_; i++) {
                flux_1g_(i) =
                    flux_1g_(i) / (xstr_(i) * vol_[i]) + qbar[i] * FPI;
            }
        } // OMP single

        cw.post_sweep();

    } // OMP Parallel

    return;
} // sweep1g
