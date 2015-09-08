#pragma once

#include "pugixml.hpp"

#include "sn_sweeper.hpp"
#include "core_mesh.hpp"
namespace mocc {
    class SnSweeper_CDD: public SnSweeper {
    public:
        SnSweeper_CDD( const pugi::xml_node &input, const CoreMesh &mesh );

        
    };
}
