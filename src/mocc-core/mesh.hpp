#pragma once

#include <iostream>
#include <memory>

#include "global_config.hpp"
#include "constants.hpp"

namespace mocc {
    /**
     * This defines a base Mesh type, which provides some basic information. For
     * now, the mesh is restricted to a structured grid of cells, which in the
     * case of the derived CoreMesh class are filled with Pin objects, but in
     * the case of the base Mesh type are more abstract in nature. In lieu of a
     * standalone coarse mesh (for things like CMFD and 2D/3D), the Mesh itself
     * provides method for interacting with homogenous regions and their
     * interface surfaces.
    */
    class Mesh {
    public:
        Mesh() { };
        Mesh( unsigned int n_reg, unsigned int n_xsreg, unsigned int nx, 
                unsigned int ny, unsigned int nz ):
            n_reg_( n_reg ),
            n_xsreg_( n_xsreg ),
            nx_( nx ),
            ny_( ny ),
            nz_( nz )
        {
            std::cout << "generating a base mesh type on its own" << std::endl;
            this->prepare_surfaces();
            return;
        }

        unsigned int n_reg() const {
            return n_reg_;
        }
        
        unsigned int nx() const {
            return nx_;
        }

        unsigned int ny() const {
            return ny_;
        }

        unsigned int nz() const {
            return nz_;
        }

        unsigned int n_pin() const {
            return nx_*ny_*nz_;
        }
        
        /**
         * Return the coarse cell index given a pin Position. Cell indexing is
         * natural in x, y z.
        */
        unsigned int coarse_cell( Position pos ) const {
            return pos.z*nx_*ny_ + pos.y*nx_ + pos.x;
        }

        /**
         * Return a coarse surface index given a cell index and Surface
         * direction. Surface indexing is more complicated than cell indexing,
         * so listen up, dear readers...
         * 
         * Imagine that you are in the bottom plane of the mesh. Start by
         * numbering all of the bottom faces of the plane, starting in the lower
         * left, then moving right and up. You will have nx_*ny_ bottom surfaces
         * indexed, from 0 to nx_*ny_-1. Now start numbering all of the x-normal
         * faces, again starting with the bottom-leftmost, proceeding right,
         * then up.  Now you should be at nx_*ny_ + (nx_+1)*ny_ - 1. Do a
         * similar thing for the y-normal faces, except proceed up first, then
         * right. Lastly, number the surfaces above you, again starting at the
         * southwest, proceeding east, then north. This entire plane is now
         * indexed. Move up to the next plane above you and repeat the process,
         * keeping in mind that the surfaces below you already have numbers.
        */
        unsigned int coarse_surf( unsigned int i, Surface surf ) const {

            return 0;
        }

    protected:
        /**
         * This method pre-computes the surface indices for each coarse cell.
         * This must be done after we know what the overall dimensions of the
         * Mesh are since there is an empty constructor, allowing for mesh
         * dimensions to be deferred (namely by the CoreMesh), these dimensions
         * need to be determined at the end.
        */
        void prepare_surfaces();

        // Total number of FSRs in the entire geometry
        unsigned int n_reg_;
        // Total number of XS regions in the entire geometry
        unsigned int n_xsreg_;
        // Numbers of pins/planes in each dimension
        unsigned int nx_;
        unsigned int ny_;
        unsigned int nz_;
    private:
        // Vector storing densely-packed coarse mesh surface indices for each
        // cell.
        VecI coarse_surf_;
    };

    typedef std::shared_ptr<Mesh> SP_Mesh_t;
}
