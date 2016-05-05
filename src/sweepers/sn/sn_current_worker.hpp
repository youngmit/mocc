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

#include "core/coarse_data.hpp"
#include "core/constants.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"

#include "util/force_inline.hpp"

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
            MOCC_FORCE_INLINE void upwind_work( const real_t *x, const real_t *y,
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
            MOCC_FORCE_INLINE void upwind_work( const real_t *x,
                    const real_t *y, const Angle &ang, int group ) {

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
            MOCC_FORCE_INLINE void current_work( real_t psi_x, real_t psi_y,
                    real_t psi_z, size_t i, Angle &ang, int group )
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
            MOCC_FORCE_INLINE void current_work( real_t psi_x, real_t psi_y,
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
            MOCC_FORCE_INLINE void set_octant( Angle ang ){
                // Configure the upwind/downwind directions based on the angle
                downwind_x_ = (ang.ox > 0.0) ? Surface::EAST : Surface::WEST;
                upwind_x_   = (ang.ox > 0.0) ? Surface::WEST : Surface::EAST;

                downwind_y_ = (ang.oy > 0.0) ? Surface::NORTH : Surface::SOUTH;
                upwind_y_   = (ang.oy > 0.0) ? Surface::SOUTH : Surface::NORTH;

                downwind_z_ = (ang.oz > 0.0) ? Surface::TOP : Surface::BOTTOM;
                upwind_z_   = (ang.oz > 0.0) ? Surface::BOTTOM : Surface::TOP;

                // Figure out which partial current to contribute
                part_x_ = (ang.ox > 0.0) ? 0 : 1;
                part_y_ = (ang.oy > 0.0) ? 0 : 1;
                part_z_ = (ang.oz > 0.0) ? 0 : 1;

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

            // Index to access for partial current for the current angle
            int part_x_;
            int part_y_;
            int part_z_;
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

            MOCC_FORCE_INLINE void upwind_work( const real_t *x,
                    const real_t *y, const real_t *z, const Angle &ang,
                    int group )
            {
                return;
            }

            MOCC_FORCE_INLINE void upwind_work( const real_t *x,
                    const real_t *y, const Angle &ang, int group )
            {
                return;
            }

            MOCC_FORCE_INLINE void current_work( real_t psi_x, real_t psi_y,
                    real_t psi_z, size_t i, const Angle &ang, int group )
            {
                return;
            }

            MOCC_FORCE_INLINE void current_work( real_t psi_x, real_t psi_y,
                    size_t i, const Angle &ang, int group )
            {
                return;
            }

            MOCC_FORCE_INLINE void set_octant( Angle ang ) {
                return;
            }
        };
    }
}

