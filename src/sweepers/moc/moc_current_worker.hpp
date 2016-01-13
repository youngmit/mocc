#pragma once

#include <cmath>

#include "angle.hpp"
#include "constants.hpp"
#include "global_config.hpp"
#include "ray.hpp"

namespace mocc {
    namespace moc {
        /**
         * This is a pure virtual class from which all valid current workers
         * should derive. There really shouldn't be any pointers of this class
         * floating around, since we never really intend to have to resolve
         * virtual method calls (use templates instead), but this provides the
         * requirements of all variants of the current worker to be used as a
         * template parameter to \ref MoCSweeper::sweep1g().
         *
         * In short, this class is provided as a means of checking that proposed
         * workers implement all of the required interfaces, and if not provide
         * a more reasonable error from the compiler than unholy template-based
         * complaints.
         */
        class CurrentWorker {
        public:
            virtual void post_ray( const ArrayF &psi5, const ArrayF &psi2,
                    const ArrayF &e_tau, const Ray &ray, int first_reg,
                    int group) = 0;

            virtual void set_angle( Angle ang, real_t spacing ) = 0;

            virtual void post_angle( int iang, int igroup ) = 0;

            virtual void set_plane( int iplane ) = 0;

            virtual void post_sweep( int igroup ) = 0;

            virtual void post_plane( int igroup ) = 0;
        };

        /**
         * This can be used as a template parameter to the \ref
         * MoCSweeper::sweep1g() method. Using this class in such a way avoids
         * the extra work needed to compute currents, and with any optimization
         * enabled, should yield code identical to a hand-written MoC sweep
         * without the current work.
         */
        class NoCurrent: public CurrentWorker {
        public:
            NoCurrent() { }
            NoCurrent( CoarseData *data, const Mesh *mesh ) { }

            /**
             * Defines work to be done following the sweep of a single ray. This
             * is useful for when you need to do something with the angular
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

            inline void set_plane( int iplane ) {
                return;
            }

            inline void post_sweep( int igroup ) {
                return;
            }

            inline void post_plane( int igroup ) {
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
        class Current: public CurrentWorker {
        public:
            Current():
                coarse_data_( nullptr ),
                mesh_( nullptr ) { }

            Current( CoarseData *data, const Mesh *mesh ):
                coarse_data_( data ),
                mesh_( mesh )
            {
                return;
            }

            inline void post_angle( int iang, int igroup ) {
                return;
            };

            inline void post_plane( int igroup ) {
                return;
            }

            inline void set_plane( int plane ) {
                plane_ = plane;
                cell_offset_ = mesh_->coarse_cell_offset( plane );
                surf_offset_ = mesh_->coarse_surf_offset( plane );
            }

            inline void set_angle( Angle ang, real_t spacing ) {
#pragma omp single
                {
                    // Scale the angle weight to sum to 4*PI
                    real_t w = ang.weight * PI;
                    // Multiply by dz so that we conform to the actual coarse
                    // mesh area.
                    real_t dz = mesh_->dz(plane_);

                    current_weights_[0] = w * ang.ox *
                        spacing/std::abs( std::cos(ang.alpha) ) *  dz;
                    current_weights_[1] = w * ang.oy *
                        spacing/std::abs( std::sin(ang.alpha) ) *  dz;
                    flux_weights_[0] = w *
                        spacing/std::abs( std::cos(ang.alpha) ) *  dz;
                    flux_weights_[1] = w *
                        spacing/std::abs( std::sin(ang.alpha) ) *  dz;
                }
#pragma omp barrier
            }

            inline void post_ray( const ArrayF &psi1, const ArrayF &psi2,
                    const ArrayF &e_tau, const Ray &ray, int first_reg,
                    int group ) {
#pragma omp critical
                {
                    ArrayB1 current =
                        coarse_data_->current(blitz::Range::all(), group);
                    ArrayB1 surface_flux =
                        coarse_data_->surface_flux(blitz::Range::all(), group);
                    size_t cell_fw = ray.cm_cell_fw()+cell_offset_;
                    size_t cell_bw = ray.cm_cell_bw()+cell_offset_;
                    int surf_fw = ray.cm_surf_fw()+surf_offset_;
                    int surf_bw = ray.cm_surf_bw()+surf_offset_;
                    size_t iseg_fw = 0;
                    size_t iseg_bw = ray.nseg();

                    Normal norm_fw = mesh_->surface_normal( surf_fw );
                    Normal norm_bw = mesh_->surface_normal( surf_bw );
                    current( surf_fw ) += psi1[iseg_fw] *
                        current_weights_[(int)norm_fw];
                    current( surf_bw ) -= psi2[iseg_bw] *
                        current_weights_[(int)norm_bw];
                    surface_flux( surf_fw ) += psi1[iseg_fw] *
                        flux_weights_[(int)norm_fw];
                    surface_flux( surf_bw ) += psi2[iseg_bw] *
                        flux_weights_[(int)norm_bw];

                    auto begin = ray.cm_data().cbegin();
                    auto end = ray.cm_data().cend();
                    for( auto crd = begin; crd != end; ++crd ) {
                        // Hopefully branch prediction saves me here.
                        if( crd->fw != Surface::INVALID ) {
                            iseg_fw += crd->nseg_fw;
                            norm_fw = surface_to_normal( crd->fw );
                            surf_fw = mesh_->coarse_surf( cell_fw, crd->fw );
                            current( surf_fw ) += psi1[iseg_fw] *
                                current_weights_[(int)norm_fw];
                            surface_flux( surf_fw ) += psi1[iseg_fw] *
                                flux_weights_[(int)norm_fw];
                        }

                        if( crd->bw != Surface::INVALID ) {
                            iseg_bw -= crd->nseg_bw;
                            norm_bw = surface_to_normal( crd->bw );
                            surf_bw = mesh_->coarse_surf( cell_bw, crd->bw );
                            current( surf_bw ) -= psi2[iseg_bw] *
                                current_weights_[(int)norm_bw];
                            surface_flux( surf_bw ) += psi2[iseg_bw] *
                                flux_weights_[(int)norm_bw];
                        }

                        cell_fw = mesh_->coarse_neighbor( cell_fw, (crd)->fw );
                        cell_bw = mesh_->coarse_neighbor( cell_bw, (crd)->bw );
                    }
                } // OMP Critical
                return;
            }


            inline void post_sweep( int igroup ) {
#pragma omp single
            {
                ArrayB1 current =
                    coarse_data_->current( blitz::Range::all(), igroup );
                ArrayB1 surface_flux =
                    coarse_data_->surface_flux( blitz::Range::all(), igroup );
                // Normalize the surface currents
                for( size_t plane=0; plane<mesh_->nz(); plane++ ) {
                    for( auto surf=mesh_->plane_surf_xy_begin(plane);
                            surf!=mesh_->plane_surf_end(plane);
                            ++surf )
                    {
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
            real_t current_weights_[2];
            real_t flux_weights_[2];

            int plane_;
            int cell_offset_;
            int surf_offset_;
        };


    }
}
