#include "mesh.hpp"

namespace mocc {
    void Mesh::prepare_surfaces() {
        // number of x-normal surfaces in each y row
        int nxsurf = nx_+1;
        // number of y-normal surfaces in each x column
        int nysurf = ny_+1;
        // number of x- and y-normal surfaces in each plane
        int nxysurf = (nx_ + 1)*ny_ + (ny_ + 1)*nx_;
        int i = 0;
        coarse_surf_ = VecI(6*nx_*ny_*nz_, 0);
        int surf_offset = 0;
        for( unsigned int iz=0; iz<nz_; iz++ ) {
            for( unsigned int iy=0; iy<ny_; iy++ ) {
                for( unsigned int ix=0; ix<nx_; ix++ ) {
                    Surface surf;
                    int cell_offset = i*6;


                    surf = Surface::EAST;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*iy + ix + 1 + nx_*ny_;

                    surf = Surface::WEST;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*iy + ix + nx_*ny_;

                    surf = Surface::NORTH;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*ny_ + nysurf*ix + iy + 1 + nx_*ny_;
                    
                    surf = Surface::SOUTH;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*ny_ + nysurf*ix + iy + nx_*ny_;

                    surf = Surface::BOTTOM;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nx_*iy + ix;

                    surf = Surface::TOP;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nx_*ny_ + nxysurf + nx_*iy + ix;

                    i++;
                }
            }
            surf_offset += nxysurf + nx_*ny_;
        }

        return;
    }

    unsigned int Mesh::coarse_surf_point( Point2 p, int cell, int (&s)[2] ) 
            const 
    {
        bool on_x = false;
        bool on_y = false;

        int ix = 0;
        for( auto &xi: x_vec_ ) {
            if( fp_equiv_abs(p.x, xi) ) {
                on_x = true;
                break;
            }
        }

        int iy = 0;
        for( auto &yi: y_vec_ ) {
            if( fp_equiv_abs(p.y, yi) ) {
                on_y = true;
                break;
            }
        }

        if( on_x && !on_y ) {
            // bottom faces (nx_*ny_) PLUS
            // previous rows (nx_+1)*iy PLUS
            // x position (ix)
            // simplified to reduce number of ops
            s[0] = nx_*(ny_+iy) + iy + ix;
            return 1;
        }

        if( on_y && !on_x ) {
            // bottom faces (nx_*ny_) PLUS
            // all x-normal surfaces (nx_+1)*ny_ PLUS
            // y position (iy)
            // simplified to reduce number of ops, fewer imuls
            s[0] = nx_*ny_ + nx_*ny_ + ny_ + iy;
            return 1;
        }

        // If we are on an x-normal and a y-normal face simultaneously, this
        // can be a potential issue for coarse ray data. For each cell in
        // the mesh to have balance, the ray flux leaving/entering a corner
        // must be accounted for, so we can't just say that the flux goes
        // directly into the diagonal neighbor, instead, we must say that it
        // goes into an adjacent neighbor first, then into the diagonal
        // neighbor, even though there are no actual ray segments in the
        // adjacent neighbor. Therefore we may need to return two surface
        // indices. in the event that we are on the border of the geometry, we
        // may only need to return one surface.
        if( on_x && on_y) {
            // For conservation, we need to be consistent with how the ray
            // should cross the corner. It needs to leave the current cell,
            // glance through the adjacent cell, and end up entering on the
            // consistent surface of the diagonal neighbor. We will use the
            // convention that the ray always goes into the x neighbor, grances
            // through, then moves in y to the diagonal neighbor.

            // Find corner of the current cell that we are looking at
            Position cellpos = this->coarse_position( cell );
            Surface corner_x = Surface::INVALID;
            Surface corner_y = Surface::INVALID;
            if( ix == cellpos.x) {
                corner_x = Surface::WEST;
            } else if ( ix == cellpos.x + 1 ) {
                corner_x = Surface::EAST;
            }
            if( iy == cellpos.y ) {
                corner_y = Surface::SOUTH;
            } else if( iy == cellpos.y + 1 ) {
                corner_y = Surface::NORTH;
            }

            assert(corner_x != Surface::INVALID);
            assert(corner_y != Surface::INVALID);

            /// So to break down the rules
            /// - on the domain boundary, only return the surface normal to the
            /// boundary. This may need to be re-addressed when we start
            /// thinking about spatial decomposition.
            /// - on the interior, go x normal first, then y normal.
            
            // Handle the boundary case
            if( (ix == 0)     && (corner_x == Surface::WEST) ) {
                s[0] = nx_*(ny_+iy) + iy + ix;
                return 1;
            }
            if( (ix == nx_+1) && (corner_x == Surface::EAST) ) {
                s[0] = nx_*(ny_+iy) + iy + ix;
                return 1;
            }
            if( (iy == 0)     && (corner_y == Surface::SOUTH) ) {
                s[0] = nx_*ny_ + nx_*ny_ + ny_ + iy;
                return 1;
            }
            if( (iy== ny_+1) && (corner_y == Surface::NORTH) ) {
                s[0] = nx_*ny_ + nx_*ny_ + ny_ + iy;
                return 1;
            }

            // So we aren't on a boundary. Handle the interior corner case
            s[0] = nx_*(ny_+iy) + iy + ix;
            s[1] = nx_*ny_ + nx_*ny_ + ny_ + iy;
            return 2;
        }

        return 0;
    }
}
