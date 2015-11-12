#pragma once

#include <vector>

#include "geom.hpp"
#include "global_config.hpp"
#include "lattice.hpp"

namespace mocc {
    class Plane {
    public:
        Plane(const std::vector<const Lattice*> &lattices, size_t nx, 
                size_t ny);
        
        const Lattice& at(size_t ix, size_t iy) const {
            return *(lattices_[ix + nx_*iy]);
        }

        /**
         * \brief Return a const pointer to the \ref PinMesh that occupies the
         * passed \ref Point2.
         *
         * \param[in] p The point at which to find a \ref PinMesh
         * \param[inout] first_reg the first region index of the returned \ref
         * Pin. The value passed in is incremented by the region offset within
         * the Plane.
         */
        const PinMesh* get_pinmesh( Point2 &p, int &first_reg) const;        

        /**
         * \brief Return a const pointer to the \ref PinMesh that is at the
         * passed \ref Position.
         */
        const PinMesh* get_pinmesh( Position pos ) const;
    
        /**
         * Return the number of solution mesh regions in the \ref Plane
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * Return the number of XS Mesh regions in the Plane
         */
        size_t n_xsreg() const {
            return n_xsreg_;
        }

        /**
         * \brief Return a vector containing the FSR volumes
         */
        VecF vols() const {
            VecF vols;
            for( auto &lat: lattices_ ) {
                for( auto &pin: *lat ) {
                    vols.insert(vols.end(), pin->vols().begin(), 
                            pin->vols().end());
                }
            }

            return vols;
        }

        /**
         * \breif Return the position of a pin, given its index.
         */
        Position pin_position( size_t ipin ) const;

    private:
        /**
         * Number of lattices in the x direction
         */
        size_t nx_;
        /**
         * Number of lattices in the y direction
         */
        size_t ny_;

        size_t n_reg_;
        size_t n_xsreg_;
        
        /**
         * Locations of \ref Lattice interfaces
         */
        VecF hx_;
        VecF hy_;

        /**
         * Local list of \ref Lattice pointers
         */
        std::vector<const Lattice*> lattices_;

        /**
         * List of the starting FSR index for each \ref Lattice in the plane
         */
        VecI first_reg_lattice_;
    };
}
