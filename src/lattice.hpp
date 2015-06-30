#pragma once

#include <memory>
#include <vector>
#include <map>

#include "pugixml.hpp"

#include "pin.hpp"
#include "global_config.hpp"
#include "arrays.hpp"

namespace mocc {
    class Lattice {
    public:
        Lattice( const pugi::xml_node &input, 
                 const std::map<int, UP_Pin_t> &pins );
        
        int id() {
            return id_;
        }
    private:
        int id_;
        int nx_;
        int ny_;
        
        // Array of pins in the lattice
        std::vector<Pin*> pins_;
        
        // Lower-left origin arrangement of pins.
        Array2D<Pin*> pins_2d_;
    };

    typedef std::shared_ptr<Lattice> SP_Lattice_t;
}
