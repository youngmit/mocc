#pragma once

#include <memory>
#include <vector>
#include <map>
#include <iostream>

#include "pugixml.hpp"

#include "pin.hpp"
#include "geom.hpp"
#include "pin_mesh_base.hpp"
#include "global_config.hpp"
#include "arrays.hpp"

namespace mocc {
    class Lattice {
    public:
        Lattice( const pugi::xml_node &input, 
                 const std::map<int, UP_Pin_t> &pins );
        
        unsigned int id() const {
            return id_;
        }

        // Number of pins in the x direction
        unsigned int nx() const {
            return nx_;
        }
        
        // Number of pins in the y direction
        unsigned int ny() const {
            return ny_;
        }

        // Total number of pins in the lattice
        unsigned int n_pin() const {
            return pins_.size();
        }

        // Return the size of the lattice along the x dimension
        real_t hx() const {
            return hx_;
        }

        // Return the size of the lattice along the y dimension
        real_t hy() const {
            return hy_;
        }

        const Pin& at( unsigned int x, unsigned int y ) const {
            assert( (0 <= x) & (x < nx_) );
            assert( (0 <= y) & (y < ny_) );
            return *pins_[y*nx_ + x];
        }

        // Return an iterator to the first Pin* in the lattice
        std::vector<Pin*>::const_iterator begin() const {
            return pins_.cbegin();
        }

        // Return an iterator past the last Pin* in the lattice
        std::vector<Pin*>::const_iterator end() const {
            return pins_.cend();
        }

        // Return a const reference to the underlying hx_vec_ array
        const VecF& hx_vec() const {
            return hx_vec_;
        }

        // Return a const reference to the underlying hy_vec_ array
        const VecF& hy_vec() const {
            return hy_vec_;
        }

        // Return the total number of regions in the lattice
        unsigned int n_reg() const {
            return n_reg_;
        }

        // Return the total number of XS regions in the lattice
        unsigned int n_xsreg() const {
            return n_xsreg_;
        }

        // Return a const reference to the Pin Mesh object located at the
        // provided point and increment the passed in first_reg by the pin's
        // first-region offset. These calls are chained from the CoreMesh ->
        // Plane -> Lattice, with each level in the geometrical heirarchy moving
        // the point to the appropriate local coordinates and offsetting the
        // first_reg value.
        const PinMesh* get_pinmesh( Point2 &p, int &first_reg ) const;

    private:
        unsigned int id_;
        unsigned int nx_;
        unsigned int ny_;
        unsigned int n_reg_;
        unsigned int n_xsreg_;
        real_t hx_;
        real_t hy_;
        VecF hx_vec_;
        VecF hy_vec_;
        VecF x_vec_;
        VecF y_vec_;
        
        // Array of pins in the lattice
        std::vector<Pin*> pins_;

        // Array of starting FSR indices for each pin in the lattice
        VecI first_reg_pin_;
    };

    typedef std::shared_ptr<Lattice> SP_Lattice_t;
}
