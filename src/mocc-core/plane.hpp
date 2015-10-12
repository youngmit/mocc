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

        const PinMesh* get_pinmesh( Point2 &p, int &first_reg) const;        
    
        /**
         * Return the number of solution mesh regions in the Plane
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

        const VecF vols() const {
            VecF vols;
            for( auto &lat: lattices_ ) {
                for( auto &pin: *lat ) {
                    vols.insert(vols.end(), pin->vols().begin(), 
                            pin->vols().end());
                }
            }

            return vols;
        }

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
        
        // Locations of lattice interfaces
        VecF hx_;
        VecF hy_;

        // Local list of lattices
        std::vector<const Lattice*> lattices_;

        // List of the starting FSR index for each lattice in the plane
        VecI first_reg_lattice_;
    };
}
