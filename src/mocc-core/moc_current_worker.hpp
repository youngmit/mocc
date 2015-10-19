#pragma once

#include "angle.hpp"
#include "constants.hpp"
#include "global_config.hpp"
#include "ray.hpp"

namespace mocc {
    namespace moc {
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
            Current():
                data_( nullptr ),
                mesh_( nullptr ) { }
            Current( CoarseData *data, const Mesh *mesh ):
                data_( data ),
                mesh_( mesh )
            {
                return;
            }

            inline void set_angle( Angle ang, real_t spacing ) {
                current_weights_[0] = ang.weight * spacing * ang.ox;
                current_weights_[1] = ang.weight * spacing * ang.oy;
                current_weights_[2] = ang.weight * spacing * ang.oz;
            }

            inline void post_ray( const ArrayF &psi1, const ArrayF &psi2, 
                    const ArrayF &e_tau, const Ray &ray, int first_reg,
                    int group ) {
                size_t cell_fw = ray.cm_cell_fw();
                size_t cell_bw = ray.cm_cell_bw();
                size_t surf_fw = ray.cm_surf_fw();
                size_t surf_bw = ray.cm_surf_bw();
                size_t iseg_fw = 0;
                size_t iseg_bw = ray.nseg();

                Normal norm_fw = mesh_->surface_normal( surf_fw );
                Normal norm_bw = mesh_->surface_normal( surf_bw );
                data_->current( surf_fw, group ) += 
                    psi1[iseg_fw] * current_weights_[(int)norm_fw];
                data_->current( surf_bw, group ) -= 
                    psi2[iseg_bw] * current_weights_[(int)norm_bw];

                auto begin = ray.cm_data().cbegin();
                auto end = ray.cm_data().cend();
                for( auto crd = begin; crd != end; ++crd ) {
                    // Hopefully branch prediction saves me here.
                    if( crd->fw != Surface::INVALID ) {
                        iseg_fw += crd->nseg_fw;
                        norm_fw = surface_to_normal( crd->fw );
                        surf_fw = mesh_->coarse_surf( cell_fw, crd->fw );
                        data_->current( surf_fw, group ) += 
                            psi1[iseg_fw] * current_weights_[(int)norm_fw];
                    }

                    if( crd->bw != Surface::INVALID ) {
                        iseg_bw -= crd->nseg_bw;
                        norm_bw = surface_to_normal( crd->bw );
                        surf_bw = mesh_->coarse_surf( cell_bw, crd->bw );
                        data_->current( surf_bw, group ) -= 
                            psi2[iseg_bw] * current_weights_[(int)norm_bw];
                    }

                    cell_fw = mesh_->coarse_neighbor( cell_fw, (crd)->fw );
                    cell_bw = mesh_->coarse_neighbor( cell_bw, (crd)->bw );
                }
                return;
            }

            inline void post_angle( int iang, int igroup ) {
                return;
            };

        private:
            CoarseData *data_;
            const Mesh *mesh_;
            real_t current_weights_[3];
        };

        /**
         * This can be used as a template parameter to the \ref
         * MoCSweeper::sweep1g() method. Using this class in such a way avoids
         * the extra work needed to compute currents, and with any optimization
         * enabled, should yield code identical to a hand-written MoC sweep
         * without the current work.
         */
        class NoCurrent {
        public:
            NoCurrent() { }
            NoCurrent( CoarseData *data, const Mesh *mesh ) { }

            /**
             * Defines work to be done following the sweep of a single ray. This
             * is useful for when you need to do something with the angluar
             * flux.
             */
            inline void post_ray( const ArrayF &psi1, const ArrayF &psi2, 
                    const ArrayF &e_tau, const Ray &ray, int first_reg,
                    int group ) {
                return;
            }

            /**
             * Defines work to be done before sweeping rays in a given angle.
             */
            inline void set_angle( Angle ang, real_t spacing ) {
                return;
            }

            /**
             * Defines work to be done after sweeping all rays in a given angle.
             */
            inline void post_angle( int iang, int igroup ) {
                return;
            };
        };
    }
}
