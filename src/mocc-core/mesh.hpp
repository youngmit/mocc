#pragma once

#include <algorithm>
#include <iostream>
#include <memory>

#include "constants.hpp"
#include "geom.hpp"
#include "global_config.hpp"

namespace mocc {
    /**
     * \page coarseraypage Coarse Ray Tracing
     * Each ray crossing a mesh corner must deposit its information
     * on one exiting face of the current cell and one entering surface of
     * the diagonal neighbor. Consistency must be maintained between
     * coincident rays of different angle, otherwise surface quantities may
     * end up with non-sensical values. A good example is when current
     * should be zero in certain symmetric situations. If the corner
     * crossings are not handled properly, non-zero current could be
     * calculated because a ray that crosses one face in one direction is
     * not being cancelled out by its sibling ray in the direction reflected
     * across that face (for instance if the reflected ray passes instead
     * through the neighboring coarse mesh surface). This would impart an
     * artificially non-zero current on both of those faces.

     * \todo include discussion of coarse ray trace peculiarities and
     * conventions.
    */

    /**
     * This defines a base Mesh type, which provides some basic information. For
     * now, the mesh is restricted to a structured grid of cells, which in the
     * case of the derived CoreMesh class are filled with Pin objects, but in
     * the case of the base Mesh type are more abstract in nature. In lieu of a
     * standalone coarse mesh (for things like CMFD and 2D/3D), the Mesh itself
     * provides method for interacting with homogenous regions and their
     * interface surfaces.
     *
     * \todo the unit test for Mesh is not really done. Its got some stuff in
     * there for my testing, but it needs some serious work.
     *
     */
    
    class Mesh {
    public:
        Mesh() { };
        Mesh( size_t n_reg, size_t n_xsreg, 
                size_t nx, size_t ny, size_t nz,
                VecF &hx, VecF &hy ):
            n_reg_( n_reg ),
            n_xsreg_( n_xsreg ),
            nx_( nx ),
            ny_( ny ),
            nz_( nz ),
            x_vec_( hx ),
            y_vec_( hy ),
            dx_vec_( hx.size()-1 ),
            dy_vec_( hy.size()-1 )
        {
            assert( std::is_sorted( x_vec_.begin(), x_vec_.end() ) );
            assert( std::is_sorted( y_vec_.begin(), y_vec_.end() ) );

            hx_ = x_vec_.back();
            for( auto &xi: x_vec_ ) {
                lines_.push_back( Line( Point2( xi, 0.0 ),
                                        Point2( xi, hx_ ) ) );
            }
            hy_ = y_vec_.back();
            for( auto &yi: y_vec_ ) {
                lines_.push_back( Line( Point2( 0.0, yi ),
                                        Point2( hx_, yi ) ) );
            }

            for( size_t i=1; i<x_vec_.size(); i++ ) {
                dx_vec_[i] = x_vec_[i] - x_vec_[i-1];
            }
            for( size_t i=1; i<y_vec_.size(); i++ ) {
                dy_vec_[i] = y_vec_[i] - y_vec_[i-1];
            }

            assert( nx == hx.size()-1 );
            assert( ny == hy.size()-1 );
            this->prepare_surfaces();
            return;
        }

        /**
         * Return the total number of regions in the computational mesh. This is
         * not neccessarily the number of pins. For a MoC/CoreMesh this is the
         * number of flat source regions.
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * Return the number of pins along the X dimension
         */
        size_t nx() const {
            return nx_;
        }

        /**
         * Return the number of pins along the Y dimension
         */
        size_t ny() const {
            return ny_;
        }

        /**
         * Return the number of planes in the Z dimension
         */
        size_t nz() const {
            return nz_;
        }
        
        /**
        * Return the total core length along the x dimension
        */
        real_t hx() const {
            return hx_;
        }

        /**
        * Return the total core length along the y dimension
        */
        real_t hy() const {
            return hy_;
        }

        /**
        * Return the pin boundary locations along the x dimension
        */
        const VecF& pin_dx() const {
            return dx_vec_;
        }

        /**
        * Return the pin boundary locations along the y dimension
        */
        const VecF& pin_dy() const {
            return dy_vec_;
        }

        /**
         * Return the total number of pin regions in the mesh. This includes
         * plane separations. This is essentially the number of coarse mesh
         * regions.
         */
        size_t n_pin() const {
            return nx_*ny_*nz_;
        }

        /**
         * Return the number of coarse surfaces.
         */
        size_t n_surf() const {
            return (nx_+1)*ny_*nz_ + (ny_+1)*nx_*nz_ + (nz_+1)*nx_*ny_;
        }

        /**
         * Return a vector containing the x-, y-, z-dimensions of the mesh.
         */
        VecI dimensions() const {
            VecI d = { (unsigned int)nx_, 
                       (unsigned int)ny_, 
                       (unsigned int)nz_ };
            return d;
        }
        
        /**
         * Return the coarse cell index given a pin Position. Cell indexing is
         * natural in x, y z.
        */
        size_t coarse_cell( Position pos ) const {
            return pos.z*nx_*ny_ + pos.y*nx_ + pos.x;
        }

        /**
         * \brief Return the coarse cell index corresponding to the Point2
         * passed.
         */
        int coarse_cell_point( Point2 p ) const;

        /**
         * Return the Position of a coarse mesh cell index.
        */
        Position coarse_position( size_t cell ) const {
            return Position( 
                    cell % nx_, 
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
        size_t coarse_surf( size_t i, Surface surf ) const {
            assert( i >= 0 );
            assert( i < this->n_pin() );
            return coarse_surf_[i*6+(int)surf];
        }

        /**
         * Determine the surface of a coarse cell that a point is on. If the
         * point is on a corner, things get weird. Read on...
         *
         * This is useful for performing the coarse ray trace, and has to
         * consider a number of concerns related to current and surface flux
         * calculations, mostly when a point is in the vicinity of a mesh
         * corner. In such a case, potentially two surfaces are crossed and both
         * must be returned. The conventions for handling this are described
         * in \ref coarseraypage.
         */
        size_t coarse_norm_point( Point2 p, int octant, Surface (&s)[2] ) const;

        /**
         * Return the cell index that a point on the boundary of the mesh should
         * be considered within.
         *
         * This follows the conventions described in
         */
        size_t coarse_boundary_cell( Point2 p, int octant ) const;

        /**
         * \brief Return the number of surfaces coincident with the passed
         * Point2 and determine the index(s) of the surface(s) crossed.
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
        int coarse_surf_point( Point2 p, int cell, int (&s)[2] ) 
            const;

        /**
         * \brief Return the neighboring coarse cell index of the passed cell
         * index in the passed direction.
         *
         * \todo Handle the case where the neighbor doesnt exist.
         */
        int coarse_neighbor( size_t cell, Surface surf) const {
            switch( surf ) {
                case Surface::NORTH:
                    return cell + nx_;
                case Surface::SOUTH:
                    return cell - nx_;
                case Surface::EAST:
                    return cell + 1;
                case Surface::WEST:
                    return cell - 1;
                case Surface::TOP:
                    return cell + nx_*ny_;
                case Surface::BOTTOM:
                    return cell - nx_*ny_;
                default:
                    return -1;
            }
        }

        /**
         * \brief Return the surface normal of the given surface.
         */
        Normal surface_normal( size_t surface ) const {
            // Number of surfaces per plane
            size_t nsurfz = nx_*ny_ + (nx_+1)*ny_ + (ny_+1)*nx_;

            if( surface % nsurfz < nx_*ny_ ) {
                return Normal::Z_NORM;
            }
            if( surface % nsurfz < nx_*ny_ + (nx_+1)*ny_ ) {
                return Normal::X_NORM;
            }
            return Normal::Y_NORM;
        }
        
        /** 
         * \brief Trace a ray through the coarse surfaces.
        */
        void trace( std::vector<Point2> &p ) const; 

    protected:
        /**
         * This method pre-computes the surface indices for each coarse cell.
         * This must be done after we know what the overall dimensions of the
         * Mesh are since there is an empty constructor, allowing for mesh
         * dimensions to be deferred (namely by the CoreMesh), these dimensions
         * need to be determined at the end.
        */
        void prepare_surfaces();

        /// Total number of FSRs in the entire geometry
        size_t n_reg_;
        /// Total number of XS regions in the entire geometry
        size_t n_xsreg_;
        // Numbers of pins/planes in each dimension
        size_t nx_;
        size_t ny_;
        size_t nz_;
        
        /// Total core size in the x dimension
        real_t hx_;

        /// Total core size in the y dimension
        real_t hy_;

        /// Total core size in the z dimension
        real_t hz_;

        /// List of pin boundaries in the x dimension (starts at 0.0)
        VecF x_vec_;

        /// List of pin boundaries in the y dimension (starts at 0.0)
        VecF y_vec_;

        /// Sequence of pin x pitches
        VecF dx_vec_;

        /// Sequence of pin y pitches
        VecF dy_vec_;
        
        /// Vector of Line objects, representing pin boundaries. This greatly
        /// simplifies the ray trace.
        std::vector<Line> lines_;

    private:
        /// Vector storing densely-packed coarse mesh surface indices for each
        /// cell.
        VecI coarse_surf_;
    };

    typedef std::shared_ptr<Mesh> SP_Mesh_t;
}
