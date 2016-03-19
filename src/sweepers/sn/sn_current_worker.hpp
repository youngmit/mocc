#pragma once

#include "coarse_data.hpp"
#include "constants.hpp"
#include "global_config.hpp"
#include "mesh.hpp"

namespace mocc {
    namespace sn {
        /**
         * This class is intended to be used as a template parameter for an
         * \ref SnSweeper sweep kernel (see \ref SnSweeperVariant::sweep_1g()
         * for an example).  When templated on this type, the sweeper kernel
         * will perform current calculations for the upwind boundary condition
         * and after sweeping each cell.
         *
         * \note Unlike in the MoC sweepers, these routines to not calculate
         * area*current. Make sure to remember to multiply by the surface areas
         * at the end of the last Sn sweep.
         */
        class Current {
        public:
            Current( CoarseData *data, const Mesh *mesh ):
                data_( data ),
                mesh_( mesh )
            {
                return;
            }

            /**
             * Store the upwind boundary condition as a contribution to the
             * coarse mesh current.
             */
            inline void upwind_work( const real_t *x, const real_t *y,
                    const real_t *z, const Angle &ang, int group ) {

                size_t nx = mesh_->nx();
                size_t ny = mesh_->ny();
                size_t nz = mesh_->nz();

                real_t w = ang.weight * HPI;

                real_t ox = ang.ox*w;
                real_t oy = ang.oy*w;
                real_t oz = ang.oz*w;

                // Configure the upwind directions based on the angle
                size_t ixx = 0;
                if( ang.ox < 0.0 ) {
                    assert( upwind_x_ == Surface::EAST );
                    ixx = mesh_->nx()-1;
                }

                size_t iyy = 0;
                if( ang.oy < 0.0 ) {
                    assert( upwind_y_ == Surface::NORTH );
                    iyy = mesh_->ny()-1;
                }

                size_t izz = 0;
                if( ang.oz < 0.0 ) {
                    assert( upwind_z_ == Surface::TOP );
                    izz = mesh_->nz()-1;
                }

                // X-normal
                for( size_t iz=0; iz<nz; iz++ ) {
                    for( size_t iy=0; iy<ny; iy++ ) {
                        Position pos( ixx, iy, iz );
                        size_t i = mesh_->coarse_cell( pos );
                        int surf = mesh_->coarse_surf( i, upwind_x_ );
#pragma omp atomic update
                        data_->current(surf, group) += ox*x[ny*iz + iy];
#pragma omp atomic update
                        data_->surface_flux(surf, group) += x[ny*iz + iy];
                    }
                }

                // Y-normal
                for( size_t iz=0; iz<nz; iz++ ) {
                    for( size_t ix=0; ix<nx; ix++ ) {
                        Position pos( ix, iyy, iz );
                        size_t i = mesh_->coarse_cell( pos );
                        int surf = mesh_->coarse_surf( i, upwind_y_ );
#pragma omp atomic update
                        data_->current(surf, group) += oy*y[nx*iz + ix];
#pragma omp atomic update
                        data_->surface_flux(surf, group) += y[nx*iz + ix];
                    }
                }

                // Z-normal
                for( size_t iy=0; iy<ny; iy++ ) {
                    for( size_t ix=0; ix<nx; ix++ ) {
                        Position pos( ix, iy, izz );
                        size_t i = mesh_->coarse_cell( pos );
                        int surf = mesh_->coarse_surf( i, upwind_z_ );
#pragma omp atomic update
                        data_->current(surf, group) += oz*z[ny*iy + ix];
#pragma omp atomic update
                        data_->surface_flux(surf, group) += z[ny*iy + ix];
                    }
                }

                return;
            }

            /**
             * Store the upwind boundary condition as a contribution to the
             * coarse mesh current (2-D version).
             */
            inline void upwind_work( const real_t *x, const real_t *y,
                    const Angle &ang, int group ) {

                size_t nx = mesh_->nx();
                size_t ny = mesh_->ny();

                real_t w = ang.weight * PI;

                real_t ox = ang.ox*w;
                real_t oy = ang.oy*w;

                // Configure the upwind directions based on the angle
                size_t ixx = 0;
                if( ang.ox < 0.0 ) {
                    assert( upwind_x_ == Surface::EAST );
                    ixx = mesh_->nx()-1;
                }

                size_t iyy = 0;
                if( ang.oy < 0.0 ) {
                    assert( upwind_y_ == Surface::NORTH );
                    iyy = mesh_->ny()-1;
                }

                // X-normal
                for( size_t iy=0; iy<ny; iy++ ) {
                    Position pos( ixx, iy, 0);
                    size_t i = mesh_->coarse_cell( pos );
                    int surf = mesh_->coarse_surf( i, upwind_x_ );
#pragma omp atomic update
                    data_->current(surf, group) += ox*x[iy];
#pragma omp atomic update
                    data_->surface_flux(surf, group) += x[iy];
                }

                // Y-normal
                for( size_t ix=0; ix<nx; ix++ ) {
                    Position pos( ix, iyy, 0 );
                    size_t i = mesh_->coarse_cell( pos );
                    int surf = mesh_->coarse_surf( i, upwind_y_ );
#pragma omp atomic update
                    data_->current(surf, group) += oy*y[ix];
#pragma omp atomic update
                    data_->surface_flux(surf, group) += y[ix];
                }

                return;
            } // Upwind work 2D

            /**
             * Store the downwind surface flux of a single cell as a
             * contribution to the coarse mesh current.
             */
            inline void current_work( real_t psi_x, real_t psi_y, real_t psi_z,
                    size_t i, Angle &ang, int group )
            {
                real_t w = ang.weight * HPI;

                real_t ox = ang.ox*w;
                real_t oy = ang.oy*w;
                real_t oz = ang.oz*w;

                // Watch out; we are assuming a direct mapping from the Sn mesh
                // index to the CM index.

                // X-normal
                {
                    int surf = mesh_->coarse_surf( i, downwind_x_ );
#pragma omp atomic update
                    data_->current( surf, group ) += psi_x*ox;
#pragma omp atomic update
                    data_->surface_flux( surf, group ) += psi_x;
                }

                // Y-normal
                {
                    int surf = mesh_->coarse_surf( i, downwind_y_ );
#pragma omp atomic update
                    data_->current( surf, group ) += psi_y*oy;
#pragma omp atomic update
                    data_->surface_flux( surf, group ) += psi_y;
                }

                // Z-normal
                {
                    int surf = mesh_->coarse_surf( i, downwind_z_ );
#pragma omp atomic update
                    data_->current( surf, group ) += psi_z*oz;
#pragma omp atomic update
                    data_->surface_flux( surf, group ) += psi_z;
                }

                return;
            }

            /**
             * Store the downwind surface flux of a single cell as a
             * contribution to the coarse mesh current. 2-D version.
             */
            inline void current_work( real_t psi_x, real_t psi_y,
                    size_t i, Angle &ang, int group )
            {
                real_t w = ang.weight * PI;

                real_t ox = ang.ox*w;
                real_t oy = ang.oy*w;

                // Watch out; we are assuming a direct mapping from the Sn mesh
                // index to the CM index.

                // X-normal
                {
                    int surf = mesh_->coarse_surf( i, downwind_x_ );
#pragma omp atomic update
                    data_->current( surf, group ) += psi_x*ox;
#pragma omp atomic update
                    data_->surface_flux( surf, group ) += psi_x;
                }

                // Y-normal
                {
                    int surf = mesh_->coarse_surf( i, downwind_y_ );
#pragma omp atomic update
                    data_->current( surf, group ) += psi_y*oy;
#pragma omp atomic update
                    data_->surface_flux( surf, group ) += psi_y;
                }

                return;
            } // current_work (2-D)

            /*
             * Configure the Current object to use the proper octant for looking
             * up surfaces.
             */
            inline void set_octant( int oct ){
                assert( (oct > 0) && (oct <=8) );

                // Configure the upwind directions based on the angle
                upwind_z_ = Surface::BOTTOM;
                downwind_z_ = Surface::TOP;
                if( oct > 4 ) {
                    upwind_z_ = Surface::TOP;
                    downwind_z_ = Surface::BOTTOM;
                }
                oct = ((oct-1) % 4) + 1;

                upwind_x_ = Surface::WEST;
                downwind_x_ = Surface::EAST;
                if( (oct == 2) || (oct == 3) ) {
                    upwind_x_ = Surface::EAST;
                    downwind_x_ = Surface::WEST;
                }

                upwind_y_ = Surface::SOUTH;
                downwind_y_ = Surface::NORTH;
                if( (oct == 3) || (oct == 4) ) {
                    upwind_y_ = Surface::NORTH;
                    downwind_y_ = Surface::SOUTH;
                }

                return;
            }


        private:
            CoarseData *data_;
            const Mesh *mesh_;
            Surface upwind_x_;
            Surface upwind_y_;
            Surface upwind_z_;
            Surface downwind_x_;
            Surface downwind_y_;
            Surface downwind_z_;

        };

        /**
         * This class is intended to be used as a template parameter for an
         * SnSweeper sweep kernel (see \ref SnSweeperVariant::sweep_1g() for an
         * example).  When templated on this type, the sweeper kernel will not
         * perform any extra work related to current calculations, which should
         * speed things up.
         */
        class NoCurrent {
        public:
            NoCurrent( CoarseData *data, const Mesh *mesh ) {
                return;
            }

            inline void upwind_work( const real_t *x, const real_t *y,
                    const real_t *z, const Angle &ang, int group )
            {
                return;
            }

            inline void upwind_work( const real_t *x, const real_t *y,
                    const Angle &ang, int group )
            {
                return;
            }

            inline void current_work( real_t psi_x, real_t psi_y, real_t psi_z,
                    size_t i, const Angle &ang, int group )
            {
                return;
            }

            inline void current_work( real_t psi_x, real_t psi_y,
                    size_t i, const Angle &ang, int group )
            {
                return;
            }

            inline void set_octant( int oct ) {
                return;
            }
        };
    }
}

