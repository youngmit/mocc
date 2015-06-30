#pragma once

#include <map>

#include "pugixml.hpp"

#include "lattice.hpp"
#include "global_config.hpp"

namespace mocc {
    class Assembly {
    public:
        Assembly( const pugi::xml_node &input, 
                  const std::map<int, Lattice> &lattices );
    private:
        int id_;
        int nz_;
        VecF hz_;
        std::map<int, Lattice*> lattices_;
    };
}
