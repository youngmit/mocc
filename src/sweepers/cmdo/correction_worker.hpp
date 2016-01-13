#pragma once

#include "core/angular_quadrature.hpp"
#include "core/coarse_data.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"
#include "core/xs_mesh_homogenized.hpp"

#include "sweepers/moc/moc_current_worker.hpp"
#include "sweepers/moc/ray.hpp"
#include "sweepers/moc/ray_data.hpp"

#include "sweepers/sn/correction_data.hpp"


namespace mocc {
    namespace cmdo {
        /**
         *
         * See documentation for \ref moc::NoCurrent for canonical documentation
         * for each of the methods.
         */
        class CurrentCorrections: public moc::Current {
        public:
            CurrentCorrections( CoarseData *coarse_data, const Mesh *mesh,
                    CorrectionData *corrections, const VectorX &qbar,
                    const ArrayF &xstr, const AngularQuadrature &ang_quad,
                    const XSMeshHomogenized &sn_xs_mesh,
                    const moc::RayData &rays ):
                moc::Current( coarse_data, mesh ),
                corrections_( corrections ),
                qbar_( qbar ),
                xstr_( xstr ),
                ang_quad_( ang_quad ),
                sn_xs_mesh_( sn_xs_mesh ),
                surf_sum_( mesh_->n_surf_plane()*2 ),
                vol_sum_( mesh_->n_cell_plane()*2 ),
                vol_norm_( mesh_->n_cell_plane() ),
                sigt_sum_( mesh_->n_cell_plane()*2 ),
                surf_norm_( mesh_->n_surf_plane()*2 ),
				rays_( rays )
            {
                return;
            }

            inline void post_ray( const ArrayF &psi1, const ArrayF &psi2,
                    const ArrayF &e_tau, const moc::Ray &ray, int first_reg,
                    int group ) {
#pragma omp critical
            {
                ArrayB1 current =
                        coarse_data_->current(blitz::Range::all(), group);
                ArrayB1 surface_flux =
                        coarse_data_->surface_flux(blitz::Range::all(), group);
                size_t cell_fw = ray.cm_cell_fw();
                size_t cell_bw = ray.cm_cell_bw();
                int surf_fw = ray.cm_surf_fw();
                int surf_bw = ray.cm_surf_bw();
                size_t iseg_fw = 0;
                size_t iseg_bw = ray.nseg();

                // Use an offset for the current contributions, but NOT for the
                // correction factors, which use plane-by-plane indexing
                Normal norm_fw = mesh_->surface_normal( surf_fw );
                Normal norm_bw = mesh_->surface_normal( surf_bw );
                current( surf_fw+surf_offset_, group ) +=
                    psi1[iseg_fw] * current_weights_[(int)norm_fw];
                current( surf_bw+surf_offset_, group ) -=
                    psi2[iseg_bw] * current_weights_[(int)norm_bw];
                surface_flux( surf_fw+surf_offset_, group ) +=
                    psi1[iseg_fw] * flux_weights_[(int)norm_fw];
                surface_flux( surf_bw+surf_offset_, group ) -=
                    psi2[iseg_bw] * flux_weights_[(int)norm_bw];

                surf_sum_[surf_fw*2+0] += psi1[iseg_fw];
                surf_sum_[surf_bw*2+1] += psi2[iseg_bw];
                surf_norm_[surf_fw*2+0] += 1.0;
                surf_norm_[surf_bw*2+1] += 1.0;

                auto begin = ray.cm_data().cbegin();
                auto end = ray.cm_data().cend();
                for( auto crd = begin; crd != end; ++crd ) {
                    // Hopefully branch prediction saves me here.
                    if( crd->fw != Surface::INVALID ) {
                        // Store forward volumetric stuff
                        for( int i=0; i<crd->nseg_fw; i++ ) {
                            int ireg = ray.seg_index(iseg_fw) + first_reg;
                            real_t xstr = xstr_[ireg];
                            real_t t = ang_.rsintheta * ray.seg_len(iseg_fw);
                            real_t fluxvol = t * qbar_[ireg] + e_tau[iseg_fw]*
                                (psi1[iseg_fw]-qbar_[ireg])/xstr;
                            vol_sum_[cell_fw*2+0] += fluxvol;
                            vol_norm_[cell_fw] += t;
                            sigt_sum_[cell_fw*2+0] += xstr*fluxvol;
                            iseg_fw++;
                        }
                        // Store FW surface stuff
                        norm_fw = surface_to_normal( crd->fw );
                        surf_fw = mesh_->coarse_surf( cell_fw, crd->fw );
                        current( surf_fw+surf_offset_, group ) +=
                            psi1[iseg_fw] * current_weights_[(int)norm_fw];
                        surface_flux( surf_fw+surf_offset_, group ) +=
                            psi1[iseg_fw] * flux_weights_[(int)norm_fw];
                        surf_sum_[surf_fw*2+0] += psi1[iseg_fw];
                        surf_norm_[surf_fw*2+0] += 1.0;
                    }

                    if( crd->bw != Surface::INVALID ) {
                        // Store backward volumetric stuff
                        for( int i=0; i<crd->nseg_bw; i++ ) {
                            iseg_bw--;
                            int ireg = ray.seg_index(iseg_bw) + first_reg;
                            real_t xstr = xstr_[ireg];
                            real_t t = ang_.rsintheta * ray.seg_len(iseg_bw);
                            real_t fluxvol = t * qbar_[ireg] + e_tau[iseg_bw] *
                                (psi2[iseg_bw+1]-qbar_[ireg])/xstr;
                            vol_sum_[cell_bw*2+1] += fluxvol;
                            sigt_sum_[cell_bw*2+1] += xstr*fluxvol;
                        }
                        // Store BW surface stuff
                        norm_bw = surface_to_normal( crd->bw );
                        surf_bw = mesh_->coarse_surf( cell_bw, crd->bw );
                        current( surf_bw+surf_offset_, group ) -=
                            psi2[iseg_bw] * current_weights_[(int)norm_bw];
                        surface_flux( surf_bw+surf_offset_, group ) -=
                            psi2[iseg_bw] * flux_weights_[(int)norm_bw];
                        surf_sum_[surf_bw*2+1] += psi2[iseg_bw];
                        surf_norm_[surf_bw*2+1] += 1.0;
                    }

                    cell_fw = mesh_->coarse_neighbor( cell_fw, (crd)->fw );
                    cell_bw = mesh_->coarse_neighbor( cell_bw, (crd)->bw );
                }
            }
                return;
            }

            inline void set_angle( Angle ang, real_t spacing ) {
                moc::Current::set_angle( ang, spacing );
                ang_ = ang;

                // Zero out all of the flux sum arrays
                surf_sum_ = 0.0;
                surf_norm_ = 0.0;
                vol_sum_ = 0.0;
                vol_norm_ = 0.0;
                sigt_sum_ = 0.0;

#pragma omp barrier
                return;
            }

            void post_angle( int iang, int igroup ) {
#pragma omp single
            {
                // Normalize the flux and sigt values and calculate
                // correction factors for the current angle/energy
                for( size_t i=0; i<vol_norm_.size(); i++ ) {
                    sigt_sum_[2*i+0] /= vol_sum_[2*i+0];
                    sigt_sum_[2*i+1] /= vol_sum_[2*i+1];
                    vol_sum_[2*i+0] /= vol_norm_[i];
                    vol_sum_[2*i+1] /= vol_norm_[i];
                }

                // Doing area-based normalization
                //surf_sum_ /= surf_norm_;

                calculate_corrections( iang, igroup );
            }
#pragma omp barrier
                return;
            }

        private:
            CorrectionData *corrections_;
            const VectorX &qbar_;
            const ArrayF &xstr_;
            const AngularQuadrature &ang_quad_;
            const XSMeshHomogenized &sn_xs_mesh_;

            ArrayF surf_sum_;
            ArrayF vol_sum_;
            ArrayF vol_norm_;
            ArrayF sigt_sum_;
            ArrayF surf_norm_;

            Angle ang_;

			const moc::RayData &rays_;

            /** \page surface_norm Surface Normalization
             * Surface normalization \todo discuss surface normalization
             */
            void calculate_corrections( size_t ang, size_t group );
        };
    }
}


