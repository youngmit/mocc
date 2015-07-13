#pragma once

#include <vector>
#include <iostream>

#include "lattice.hpp"
#include "geom.hpp"
#include "global_config.hpp"

namespace mocc {
    class Plane {
    public:
        Plane(const std::vector<const Lattice*> &lattices, unsigned int nx, 
                unsigned int ny);
        
        const Lattice& at(unsigned int ix, unsigned int iy) const {
            return *(lattices_[ix + nx_*iy]);
        }

        const PinMesh* get_pinmesh( Point2 &p, int &first_reg) const;        

        unsigned int n_reg() const {
            return n_reg_;
        }

        unsigned int n_xsreg() const {
            return n_xsreg_;
        }

        const VecF vol() const {
            VecF vol;
            for( auto &lat: lattices_ ) {
                for( auto &pin: *lat ) {
                    vol.insert(vol.end(), pin->vol().begin(), pin->vol().end());
                }
            }

            return vol;
        }
    private:
        // Plane dimensions in lattices
        unsigned int nx_;
        unsigned int ny_;

        unsigned int n_reg_;
        unsigned int n_xsreg_;
        
        // Locations of lattice interfaces
        VecF hx_;
        VecF hy_;
        // Local list of lattices
        std::vector<const Lattice*> lattices_;
        // List of the starting FSR index for each lattice in the plane
        VecI first_reg_lattice_;
    };
}
