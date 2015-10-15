#pragma once

#include "global_config.hpp"
#include "mesh.hpp"
#include "coarse_data.hpp"

namespace mocc {
    namespace sn {
        /**
         * This class is intended to be used as a template parameter for an
         * SnSweeper sweep kernel (see SnSweeper::sweep_dd() for an example).
         * When templated on this type, the sweeper kernel will perform current
         * calculations for the upwind boundary condition and after sweeping
         * each cell.
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
            inline void upwind_work( const ArrayF &x, const ArrayF &y, 
                    const ArrayF &z, const Angle &ang, int group ) {

                size_t nx = mesh_->nx();
                size_t ny = mesh_->ny();
                size_t nz = mesh_->nz();

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
                        size_t surf = mesh_->coarse_surf( i, upwind_x_ );
                        data_->current(surf, group) += ang.ox*x[ny*iz + iy];
                    }
                }
                
                // Y-normal
                for( size_t iz=0; iz<nz; iz++ ) {
                    for( size_t ix=0; ix<nx; ix++ ) {
                        Position pos( ix, iyy, iz );
                        size_t i = mesh_->coarse_cell( pos );
                        size_t surf = mesh_->coarse_surf( i, upwind_y_ );
                        data_->current(surf, group) += ang.oy*y[nx*iz + ix];
                    }
                }

                // Z-normal
                for( size_t iy=0; iy<ny; iy++ ) {
                    for( size_t ix=0; ix<nx; ix++ ) {
                        Position pos( ix, iy, izz );
                        size_t i = mesh_->coarse_cell( pos );
                        size_t surf = mesh_->coarse_surf( i, upwind_z_ );
                        data_->current(surf, group) += ang.oz*z[ny*iy + ix];
                    }
                }
                
                return;
            }

            /** 
             * Store the downwind surface flux of a single cell as a
             * contribution to the coarse mesh current.
             */
            inline void current_work( real_t psi_x, real_t psi_y, real_t psi_z,
                    size_t i, Angle &ang, int group ) {
                // Watch out; we are assuming a direct mapping from the Sn mesh
                // index to the CM index.

                // X-normal
                {
                    size_t surf = mesh_->coarse_surf( i, downwind_x_ );
                    data_->current( surf, group ) += psi_x*ang.ox;
                }

                // Y-normal
                {
                    size_t surf = mesh_->coarse_surf( i, downwind_y_ );
                    data_->current( surf, group ) += psi_y*ang.oy;
                }
                
                // Z-normal
                {
                    size_t surf = mesh_->coarse_surf( i, downwind_z_ );
                    data_->current( surf, group ) += psi_z*ang.oz;
                }

                return;
            }

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
         * SnSweeper sweep kernel (see SnSweeper::sweep_dd() for an example).
         * When templated on this type, the sweeper kernel will not perform any
         * extra work related to current calculations, which should speed things
         * up.
         */
        class NoCurrent {
        public:
            NoCurrent( CoarseData *data, const Mesh *mesh ) {
                return;
            }
            
            inline void upwind_work( const ArrayF &x, const ArrayF &y, 
                    const ArrayF &z, const Angle &ang, int group )
            {
                return;
            }

            inline void current_work( real_t psi_x, real_t psi_y, real_t psi_z,
                    size_t i, Angle &ang, int group )
            {
                return;
            }

            inline void set_octant( int oct ) {
                return;
            }
        };
    }
}

