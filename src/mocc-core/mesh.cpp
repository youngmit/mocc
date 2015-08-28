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


                    surf = EAST;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*iy + ix + 1 + nx_*ny_;

                    surf = WEST;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*iy + ix + nx_*ny_;

                    surf = NORTH;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*ny_ + nysurf*ix + iy + 1 + nx_*ny_;
                    
                    surf = SOUTH;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nxsurf*ny_ + nysurf*ix + iy + nx_*ny_;

                    surf = BOTTOM;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nx_*iy + ix;

                    surf = TOP;
                    coarse_surf_[cell_offset+(int)surf] = surf_offset + 
                        nx_*ny_ + nxysurf + nx_*iy + ix;

                    i++;
                }
            }
            surf_offset += nxysurf + nx_*ny_;
        }

        return;
    }
}
