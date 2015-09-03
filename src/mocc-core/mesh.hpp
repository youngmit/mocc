#pragma once

#include <iostream>
#include <memory>

#include "global_config.hpp"
#include "constants.hpp"
#include "geom.hpp"

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
        Mesh( unsigned int n_reg, unsigned int n_xsreg, 
                unsigned int nx, unsigned int ny, unsigned int nz,
                VecF &hx, VecF &hy ):
            n_reg_( n_reg ),
            n_xsreg_( n_xsreg ),
            nx_( nx ),
            ny_( ny ),
            nz_( nz ),
            x_vec_( hx ),
            y_vec_( hy )
        {
            std::cout << "generating a base mesh type on its own" << std::endl;
            hx_ = x_vec_.back();
            hy_ = y_vec_.back();
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
         * Return the Position of a coarse mesh cell index.
        */
        Position coarse_position( unsigned int cell ) const {
            return Position( 
                    cell % ny_, 
                    (cell % (nx_*ny_)) / nx_,
                    cell / (nx_*ny_) );
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
            return coarse_surf_[i*6+(int)surf];
        }

        /**
         * \brief Return the number of surfaces coincident with the passed
         * Point2.
         *
         * \param[in] p the Point2 to check for surface.
         * \param[in] cell the index of a cell that the point should be
         * bordering.
         * \param[out] s and array[2] to store the indices of the surface(s)
         * crossed.
         *
         * This will return the number of surfaces that a point lies on, either
         * 0, 1 or 2.
         *
         * Since this is primarily used for ray tracing, it is fundamentally 2-D
         * in nature, so the passed cell index will be %-ed by the number of CM
         * cells per plane, and the surface indices are returned as if they were
         * in the bottom-most plane; the code using the resulting indices is
         * therefore required to offset them to the appropriate plane.
        */
        unsigned int coarse_surf_point( Point2 p, int cell, int (&s)[2] ) 
            const;

        /**
         * \brief Return the surface normal of the given surface.
         */
        Surface surface_normal( int surface ) const {
            /// \todo implement this
            return Surface::NORTH;
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
        
        // Total core size in the x dimension
        float_t hx_;

        // Total core size in the y dimension
        float_t hy_;

        // Total core size in the z dimension
        float_t hz_;

        // List of pin boundaries in the x dimension (starts at 0.0)
        VecF x_vec_;

        // List of pin boundaries in the y dimension (starts at 0.0)
        VecF y_vec_;

    private:
        // Vector storing densely-packed coarse mesh surface indices for each
        // cell.
        VecI coarse_surf_;
    };

    typedef std::shared_ptr<Mesh> SP_Mesh_t;
}
