#pragma once

#include <map>
#include <vector>
#include <memory>

#include "pugixml.hpp"

#include "lattice.hpp"
#include "global_config.hpp"

namespace mocc {
    class Assembly {
    public:
        Assembly( const pugi::xml_node &input, 
                  const std::map<int, Lattice> &lattices );

        ~Assembly();

        int id() const {
            return id_;
        }
    private:
        int id_;
        int nz_;
        VecF hz_;
        std::vector<const Lattice*> lattices_;
    };

    typedef std::unique_ptr<Assembly> UP_Assembly_t;
}
