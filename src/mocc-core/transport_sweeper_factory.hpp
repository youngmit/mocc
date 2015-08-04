#pragma once

#include "pugixml.hpp"

#include "transport_sweeper.hpp"

namespace mocc {
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
        const CoreMesh &mesh );
}
