#pragma once

#include <map>

#include "pugixml.hpp"

#include "arrays.hpp"
#include "assembly.hpp"

namespace mocc {
    class Core {
    public:
        Core();
        Core( const pugi::xml_node &input, 
              const std::map<int, UP_Assembly_t> &assemblies);
        ~Core();

    private:
        // Core dimensions
        int nx_;
        int ny_;
        // 2D array of assemblies
        Array2D<Assembly*> assemblies_;
    };
}
