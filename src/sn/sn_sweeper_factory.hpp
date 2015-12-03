#pragma once

#include "pugixml.hpp"

#include "core_mesh.hpp"
#include "sn_sweeper.hpp"

namespace mocc {
    UP_Sweeper_t SnSweeperFactory( const pugi::xml_node &input,
            const CoreMesh &mesh );
}
