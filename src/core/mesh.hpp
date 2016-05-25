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

#include <array>
#include <algorithm>
#include <memory>

#include "constants.hpp"
#include "error.hpp"
#include "geometry/geom.hpp"
#include "global_config.hpp"
#include "position.hpp"

namespace mocc {
    /**
     * \page coarseraypage Coarse Ray Tracing
     * Each ray crossing a mesh corner must deposit its information on one
     * exiting face of the current cell and one entering surface of the diagonal
     * neighbor. Consistency must be maintained between coincident rays of
     * different angle, otherwise surface quantities may end up with
     * nonsensical values. A good example is when current should be zero in
     * certain symmetric situations. If the corner crossings are not handled
     * properly, non-zero current could be calculated because a ray that crosses
     * one face in one direction is not being cancelled out by its sibling ray
     * in the direction reflected across that face (for instance if the
     * reflected ray passes instead through the neighboring coarse mesh
     * surface). This would impart an artificially non-zero current on both of
     * those faces.

     * \todo include discussion of coarse ray trace peculiarities and
     * conventions.
    */

    /**
     * This defines a base \ref Mesh type, which provides some basic
     * information. For now, the mesh is restricted to a structured grid of
     * cells, which in the case of the derived \ref CoreMesh class are filled
     * with \ref Pin objects, but in the case of the base \ref Mesh type are
     * more abstract in nature. In lieu of a standalone coarse mesh (for things
     * like CMFD and 2D/3D), the Mesh itself provides method for interacting
     * with homogeneous regions and their interface surfaces.
     *
     * \todo the unit test for Mesh is not really done. Its got some stuff in
     * there for my testing, but it needs some serious work.
     */
    class Mesh {
    public:
        typedef std::array< std::array< Boundary, 2 >, 3 > BCArray_t;

        Mesh() { };

        /**
         * Disable the copy constructor
         */
        Mesh( const Mesh &other ) = delete;

        /**
         * \brief Construct a \ref Mesh using cell boundaries specified
         * externally.
         */
        Mesh( size_t n_reg, size_t n_xsreg,
                VecF &hx, VecF &hy, VecF &hz, std::array<Boundary, 6> bc );

        /**
         * \brief Return the total number of regions in the computational mesh.
         *
         * This is not necessarily the number of pins. For a MoC/CoreMesh this
         * is the number of flat source regions.
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * \brief Return the number of pins along the X dimension
         */
        size_t nx() const {
            return nx_;
        }

        /**
         * \brief Return the number of pins along the Y dimension
         */
        size_t ny() const {
            return ny_;
        }

        /**
         * \brief Return the number of planes in the Z dimension
         */
        size_t nz() const {
            return nz_;
        }

        /**
        * \brief Return a vector containing the core boundary conditions.
        */
        const std::array<Boundary, 6>& boundary() const {
            return bc_;
        }

        /**
         * \brief Return the boundary condition associated with the given \ref
         * Surface of the mesh
         */
        Boundary boundary_condition( Surface surf ) const {
            assert((int)surf >= (int)Surface::EAST);
            assert((int)surf <= (int)Surface::BOTTOM);

            return bc_[(int)surf];
        }

        /**
         * \brief Return a nifty array of the boundary conditions, organized by
         * normal direction and sense.
         *
         * This returns a 2-dimensional array of the boundary conditions,
         * organized in such a way that makes BC determination for a given
         * normal direction and directional sense. For example \c
         * bc[Normal::X_NORM][0] would give the boundary condition on the west
         * bounday, \c bc[Normal::Z_NORM][1] would give the boundary condition
         * on the top.
         */
        BCArray_t boundary_array() const {
            BCArray_t bc;

            bc[(int)Normal::X_NORM][0] = bc_[(int)Surface::WEST];
            bc[(int)Normal::X_NORM][1] = bc_[(int)Surface::EAST];
            bc[(int)Normal::Y_NORM][0] = bc_[(int)Surface::SOUTH];
            bc[(int)Normal::Y_NORM][1] = bc_[(int)Surface::NORTH];
            bc[(int)Normal::Z_NORM][0] = bc_[(int)Surface::BOTTOM];
            bc[(int)Normal::Z_NORM][1] = bc_[(int)Surface::TOP];

            return bc;
        }

        /**
        * \brief Return the total core length along the x dimension
        */
        real_t hx_core() const {
            return hx_;
        }

        /**
        * \brief Return the total core length along the y dimension
        */
        real_t hy_core() const {
            return hy_;
        }

        /**
         * \brief Return the pin/coarse cell thickness in the x dimension at the
         * specified x position.
         */
        inline real_t dx( size_t ix ) const {
            return dx_vec_[ix];
        }

        /**
         * \brief Return the pin/coarse cell thickness in the y dimension at the
         * specified y position.
         */
        inline real_t dy( size_t iy ) const {
            return dy_vec_[iy];
        }

        /**
         * \brief Return the pin/coarse cell thickness in the z dimension at the
         * specified z position.
         */
        inline real_t dz( size_t iz ) const {
            return dz_vec_[iz];
        }

        /**
         * \brief Return the indexed plane boundary.
         *
         * The 0-th value should be 0.0, and the last value should be to total
         * height of the domain.
         */
        inline real_t z( int iz ) const {
            return z_vec_[iz];
        }

        /**
        * \brief Return the pin boundary locations along the x dimension
        */
        const VecF& pin_dx() const {
            return dx_vec_;
        }

        /**
        * \brief Return the pin boundary locations along the y dimension
        */
        const VecF& pin_dy() const {
            return dy_vec_;
        }

        real_t coarse_volume( size_t cell ) const {
            return vol_[cell];
        }

        const VecF& coarse_volume() const {
            return vol_;
        }

        /**
         * \brief Return the cell thickness in the direction indicated.
         *
         * \param cell the index of the cell
         * \param norm the \ref Normal direction
         *
         * Given a cell index and a \ref Normal, returns the thickness of the
         * indexed cell in the specified direction. For instance, if \ref
         * Normal::X_NORM is passed, this returns the thickness of the indexed
         * cell in the x direction.
         *
         */
        real_t cell_thickness( size_t cell, Normal norm ) const {
            assert( cell < this->n_pin() );
            assert( (norm == Normal::X_NORM) || (norm == Normal::Y_NORM) ||
                    (norm == Normal::Z_NORM) );
            auto pos = this->coarse_position( cell );
            real_t h = 0;
            switch( norm ) {
            case Normal::X_NORM:
                h = dx_vec_[pos.x];
            case Normal::Y_NORM:
                h = dy_vec_[pos.y];
            case Normal::Z_NORM:
                h = dz_vec_[pos.z];
            }
            return h;
        }

        /**
         * \brief Return the total number of pin regions in the mesh.
         *
         * This includes plane separations. This is essentially the number of
         * coarse mesh regions.
         */
        size_t n_pin() const {
            return nx_*ny_*nz_;
        }

        /**
         * \brief Return the number of coarse surfaces.
         */
        size_t n_surf() const {
            return (nx_+1)*ny_*nz_ + (ny_+1)*nx_*nz_ + (nz_+1)*nx_*ny_;
        }

        /**
         * \brief Return a vector containing the x-, y-, z-dimensions of the
         * mesh.
         */
        VecI dimensions() const {
            VecI d = { (int)nx_,
                       (int)ny_,
                       (int)nz_ };
            return d;
        }

        /**
         * \brief Return the lowest coarse cell index in a given plane
         */
        int plane_cell_begin( size_t plane ) const {
            return nx_*ny_*plane;
        }

        /**
         * \brief Return the highest coarse cell index in a given plane, plus 1
         */
        int plane_cell_end( size_t plane ) const {
            return nx_*ny_*(plane+1);
        }

        /**
         * \brief Return the lowest coarse surface index in a given plane.
         *
         * It should be considered possible for the return value to be of any
         * direction normal (e.g. X, Y, or Z). Under the current indexing
         * scheme, it will always be Z-normal, but this could change.
         */
        int plane_surf_begin( size_t plane ) const {
            return n_surf_plane_*plane;
        }

        /**
         * \brief Return the highest coarse surface index in a given plane, plus
         * 1
         */
        int plane_surf_end( size_t plane ) const {
            return n_surf_plane_*(plane+1);
        }

        /**
         * \brief Return the lowest coarse surface index of the x- and y-normal
         * surfaces in a given plane.
         *
         * Iterating from here to \ref plane_surf_end() is safe, since the x/y
         * surfaces are guaranteed contiguous.
         */
        int plane_surf_xy_begin( size_t plane ) const {
            return n_surf_plane_*plane + nx_*ny_;
        }

        /**
         * \brief Return the number of coarse cells per plane
         */
        size_t n_cell_plane() const {
            return nx_*ny_;
        }

        /**
         * \brief Return the number of surfaces in each plane
         *
         * The "number of surfaces in each plane" can be a little confusing.
         * This really means the number of surfaces you need to offset one
         * surface index by to end up with the same plane-local surface, but in
         * the plane above. This means that \c n_surf_plane does not include the
         * "top" surfaces of each plane. For a concrete example, a 3x3 array of
         * regions per plane would yield an \c n_surf_plane of 33.
         */
        size_t n_surf_plane() const {
            return n_surf_plane_;
        }

        /**
         * \brief Return the coarse cell index given a pin \ref Position.
         *
         * Cell indexing is natural in x, y z.
        */
        inline int coarse_cell( Position pos ) const {
            return pos.z*nx_*ny_ + pos.y*nx_ + pos.x;
        }

        /**
         * \brief Return the coarse cell index corresponding to the \ref Point2
         * passed.
         */
        int coarse_cell_point( Point2 p ) const;

        /**
         * \brief Return the \ref Position of a coarse mesh cell index.
        */
        Position coarse_position( size_t cell ) const {
            return Position(
                    cell % nx_,
                    (cell % (nx_*ny_)) / nx_,
                    cell / (nx_*ny_) );
        }

        /**
         * \brief Return a coarse surface index given a cell index and Surface
         * direction.
         *
         * Surface indexing is more complicated than cell indexing, so listen
         * up, dear readers...
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
        int coarse_surf( size_t i, Surface surf ) const {
            assert( i >= 0 );
            assert( i < this->n_pin() );
            return coarse_surf_[i*6+(int)surf];
        }

        /**
         * \brief Determine the surface of a coarse cell that a point is on.
         *
         * If the point is on a corner, things get weird. Read on...
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
         * \brief Return the cell index that a point on the boundary of the mesh
         * should be considered within.
         *
         * This follows the conventions described in \ref coarseraypage.
         */
        int coarse_boundary_cell( Point2 p, int octant ) const;

        /**
         * \brief Return the number of surfaces coincident with the passed
         * \ref Point2 and determine the index(s) of the surface(s) crossed.
         *
         * \param[in] p the \ref Point2 to check for surface.
         * \param[in] cell the index of a cell that the point should be
         * bordering.
         * \param[out] s an std::array<int, 2> to store the indices of the
         * surface(s) crossed.
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
        int coarse_surf_point( Point2 p, int cell,
                std::array<int, 2> &s ) const;

        /**
         * \brief Return the coarse cell indices straddling the surface
         * indicated.
         *
         * \param surf the index of the coarse surface to find neighboring
         * cells.
         *
         * \note The cells are always given in the order of increasing
         * position. This makes it easy to tell which cell is to the "right"
         * of the other. The first cell in the pair is always "left" of the
         * surface and the second is always "right." Positive current always
         * flows "right."
         *
         */
        std::pair<int, int> coarse_neigh_cells( size_t surf ) const {
            std::pair<int, int> cells;

            switch( this->surface_normal(surf) ) {
            case Normal::X_NORM:
                {
                    // iz and iy are the cell positions, ix is the surface
                    // position.
                    int iz = surf / n_surf_plane_;
                    surf -= iz*n_surf_plane_;
                    int iy = (surf - nx_*ny_) / (nx_+1);
                    int ix = (surf - nx_*ny_) % (nx_+1);

                    if( ix > 0 ) {
                        cells.first = this->coarse_cell(Position(ix-1, iy, iz));
                    } else {
                        cells.first = -1;
                    }
                    if( ix < nx_ ) {
                        cells.second = this->coarse_cell(Position(ix, iy, iz));
                    } else {
                        cells.second = -1;
                    }
                }
                break;
            case Normal::Y_NORM:
                {
                    // iz and ix are the cell positions, iy is the surface
                    // position.
                    int iz = surf / n_surf_plane_;
                    surf -= iz*n_surf_plane_;
                    surf -= nx_*ny_ + (nx_+1)*ny_;
                    int iy = surf % (ny_+1);
                    int ix = surf / (ny_+1);

                    if( iy > 0 ) {
                        cells.first = this->coarse_cell(Position(ix, iy-1, iz));
                    } else {
                        cells.first = -1;
                    }
                    if( iy < ny_ ) {
                        cells.second = this->coarse_cell(Position(ix, iy, iz));
                    } else {
                        cells.second = -1;
                    }
                }
                break;
            case Normal::Z_NORM:
                {
                    // ix and iy are the cell positions, iz is the surface
                    // position.
                    int iz = surf / n_surf_plane_;
                    surf -= iz*n_surf_plane_;
                    int iy = surf / nx_;
                    int ix = surf % nx_;

                    if( iz > 0 ) {
                        cells.first = this->coarse_cell(Position(ix, iy, iz-1));
                    } else {
                        cells.first = -1;
                    }
                    if( iz < nz_ ) {
                        cells.second = this->coarse_cell(Position(ix, iy, iz));
                    } else {
                        cells.second = -1;
                    }
                }
                break;
            }
        return cells;
        }

        /**
         * \brief Return an index offset to the zero-th coarse cell in a given
         * plane.
         */
        int coarse_cell_offset( size_t plane ) const {
            return nx_*ny_*plane;
        }

        /**
         * \brief Return an index offset to the zero-th coarse surface in a
         * given plane.
         *
         * Plane indexing is wonkier than cell indexing, espectially when
         * talking about planes, since we need to take into consideration the
         * surfaces between the planes. Since surface indexing starts with the
         * bottom-facing Z-normal faces, surface "zero" for any plane is a
         * surface looking "down" into the plane below. Stuff that wants to
         * interact with X- and Y-normal faces (e.g. 2-D MoC) should expect
         * this.
         */
        int coarse_surf_offset( size_t plane ) const {
            return n_surf_plane_*plane;
        }

        /**
         * \brief Return the surface index between two cells and the \ref
         * Surface of \p cell1 that makes the interface.
         *
         * This returns a \c std::pair containing the surface index and a \ref
         * Surface enumeration. The surface index is self-explainatory. The \ref
         * Surface enumeration is a little more complicated; it is the \ref
         * Surface of \p cell1 that makes the interface between \p cell1 and \p
         * cell2. If \p cell2 is to the north of \p cell1, then the \ref Surface
         * will be \ref Surface::NORTH. If the cell indices were passed in the
         * reverse order, the \ref Surface would be \ref Surface::SOUTH.
         * Therefore, it is important which order the cells are provided.
         */
        std::pair< size_t, Surface > coarse_interface( size_t cell1,
                size_t cell2 ) const {
            // Stupid search
            for( auto is: AllSurfaces ) {
                if( this->coarse_neighbor( cell1, is ) == (int)cell2 ) {
                    return std::pair< size_t, Surface >(
                            this->coarse_surf( cell1, is ), is );
                }
            }

            throw EXCEPT("Cells do not appear to be neighbors");
        }
        /**
         * \brief Return the surface area of the indexed cell and direction.
         */
        real_t coarse_area( size_t cell, Surface surf ) const {
            auto pos = this->coarse_position( cell );
            real_t area = 0.0;
            switch( surf ) {
                case( Surface::EAST ):
                    area = dy_vec_[pos.y] * dz_vec_[pos.z];
                    break;
                case( Surface::WEST ):
                    area = dy_vec_[pos.y] * dz_vec_[pos.z];
                    break;
                case( Surface::NORTH ):
                    area = dx_vec_[pos.x] * dz_vec_[pos.z];
                    break;
                case( Surface::SOUTH ):
                    area = dx_vec_[pos.x] * dz_vec_[pos.z];
                    break;
                case( Surface::TOP ):
                    area = dx_vec_[pos.x] * dy_vec_[pos.y];
                    break;
                case( Surface::BOTTOM ):
                    area = dx_vec_[pos.x] * dy_vec_[pos.y];
                    break;
                default:
                    assert(false);
                    area = -1;
            }
            return area;
        }

        /**
         * \brief Return the surface area of the indexed coarse surface.
         */
        real_t coarse_area( size_t surf ) const {
            real_t area;
            switch( this->surface_normal(surf) ) {
            case Normal::X_NORM:
                {
                    // iz and iy are the cell positions, ix is the surface
                    // position.
                    size_t iz = surf / n_surf_plane_;
                    surf -= iz*n_surf_plane_;
                    size_t iy = (surf - nx_*ny_) / (nx_+1);

                    area = dy_vec_[iy] * dz_vec_[iz];
                }
                break;
            case Normal::Y_NORM:
                {
                    // iz and ix are the cell positions, iy is the surface
                    // position.
                    size_t iz = surf / n_surf_plane_;
                    surf -= iz*n_surf_plane_;
                    surf -= nx_*ny_ + (nx_+1)*ny_;
                    size_t ix = surf / (ny_+1);

                    area = dx_vec_[ix] * dz_vec_[iz];
                }
                break;
            case Normal::Z_NORM:
                {
                    // ix and iy are the cell positions, iz is the surface
                    // position.
                    size_t iz = surf / n_surf_plane_;
                    surf -= iz*n_surf_plane_;
                    size_t iy = surf / nx_;
                    size_t ix = surf % nx_;

                    area = dx_vec_[ix] * dy_vec_[iy];
                }
                break;
            }
            return area;
        }

        /**
         * \brief Return the neighboring coarse cell index of the passed cell
         * index in the passed direction.
         *
         * If the neighbor is out of bounds, return -1;
         */
        int coarse_neighbor( size_t cell, Surface surf) const {
            auto pos = this->coarse_position( cell );
            switch( surf ) {
                case Surface::NORTH:
                    if( pos.y < ny_-1 ) {
                        return cell + nx_;
                    } else {
                        return -1;
                    }
                case Surface::SOUTH:
                    if( pos.y > 0 ) {
                        return cell - nx_;
                    } else {
                        return -1;
                    }
                case Surface::EAST:
                    if( pos.x < nx_-1 ) {
                        return cell + 1;
                    } else {
                        return -1;
                    }
                case Surface::WEST:
                    if( pos.x > 0 ) {
                        return cell - 1;
                    } else {
                        return -1;
                    }
                case Surface::TOP:
                    if( pos.z < nz_-1 ) {
                        return cell + nx_*ny_;
                    } else {
                        return -1;
                    }
                case Surface::BOTTOM:
                    if( pos.z > 0 ) {
                        return cell - nx_*ny_;
                    } else {
                        return -1;
                    }
                default:
                    return -5;
            }
        }

        /**
         * \brief Return the plane index for the given axial position
         *
         * \param z the axial location
         * \param oz a direction in the z direction to be used to disabmiguate
         * between planes, if the value \p z lies directly on a plane interface.
         * Defaults to zero, and gives the following behavior:
         *  - If zero: Throw an exception if \p z lies on a plane interface.
         *  - If negative: Return the lower plane index if \p z lies on a plane
         *  interface.
         *  - If positive: Return the higher plane index if \p z lies on a plane
         *  interface.
         */
        int plane_index( real_t z, real_t oz=0.0 ) const {
            int iz = std::distance(z_vec_.begin(),
                    std::lower_bound(z_vec_.begin(), z_vec_.end(), z,
                                     fuzzy_lt));

            if( fp_equiv_abs(z, z_vec_[iz])) {
                if( oz == 0.0 ) {
                    throw EXCEPT("Ambiguous plane index, without valid "
                            "z-direction.");
                } else {
                    iz = (oz > 0.0) ? iz+1 : iz;
                }
            }

            return iz-1;
        }

        /**
         * \brief Return the surface normal of the given surface.
         */
        Normal surface_normal( int surface ) const {
            // Number of surfaces per plane
            int nsurfz = nx_*ny_ + (nx_+1)*ny_ + (ny_+1)*nx_;

            if( (surface % nsurfz) < nx_*ny_ ) {
                return Normal::Z_NORM;
            }
            if( (surface % nsurfz) < (nx_*ny_ + (nx_+1)*ny_) ) {
                return Normal::X_NORM;
            }
            return Normal::Y_NORM;
        }

        /**
         * \brief Trace a ray through the coarse surfaces.
        */
        void trace( std::vector<Point2> &p ) const;

        /**
         * \brief Return a const reference to the collection of \ref Line
         * objects.
         */
        const std::vector<Line> lines() const {
            return lines_;
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

        /// Total number of FSRs in the entire geometry
        size_t n_reg_;
        /// Total number of XS regions in the entire geometry
        size_t n_xsreg_;
        // Numbers of pins/planes in each dimension
        int nx_;
        int ny_;
        int nz_;

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

        /// Sequence of plane boundaries in the z dimension (starts at 0.0)
        VecF z_vec_;

        /// Sequence of pin x pitches
        VecF dx_vec_;

        /// Sequence of pin y pitches
        VecF dy_vec_;

        /// Sequence of plane heights
        VecF dz_vec_;

        /// Coarse cell volumes
        VecF vol_;

        /// Vector of \ref Line objects, representing pin boundaries. This
        /// greatly simplifies the ray trace.
        std::vector<Line> lines_;

        /// Number of surfaces per plane (doesn't consider the top surface)
        size_t n_surf_plane_;

        /// Boundary condition for each side of the mesh
        std::array<Boundary, 6> bc_;

    private:
        /// Vector storing densely-packed coarse mesh surface indices for each
        /// cell.
        VecI coarse_surf_;
    };

    typedef std::shared_ptr<Mesh> SP_Mesh_t;
}
