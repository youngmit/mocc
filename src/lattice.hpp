#pragma once

#include <memory>
#include <vector>
#include <map>
#include <iostream>

#include "pugixml.hpp"

#include "pin.hpp"
#include "global_config.hpp"
#include "arrays.hpp"

namespace mocc {
    class Lattice {
    public:
        Lattice( const pugi::xml_node &input, 
                 const std::map<int, UP_Pin_t> &pins );
        
        int id() const {
            return id_;
        }

        int nx() const {
            return nx_;
        }

        int ny() const {
            return ny_;
        }

        float_t hx() const {
            return hx_;
        }

        float_t hy() const {
            return hy_;
        }

        Pin* at( int x, int y ) const {
            assert( 0 <= x & x < nx_ );
            assert( 0 <= y & y < ny_ );
            return pins_[y*nx_ + x];
        }
    private:
        int id_;
        int nx_;
        int ny_;
        float_t hx_;
        float_t hy_;
        VecF hx_vec_;
        VecF hy_vec_;
        
        // Array of pins in the lattice
        std::vector<Pin*> pins_;
    };

    typedef std::shared_ptr<Lattice> SP_Lattice_t;
}
