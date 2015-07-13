#pragma once

#include "pugixml.hpp"

#include "transport_sweeper.hpp"
#include "moc_sweeper.hpp"
//#include "planar_synthesis_2d3d_sweeper.hpp"

namespace mocc {
    TransportSweeper* TransportSweeperFactory( const pugi::xml_node &input,
        const CoreMesh &mesh );
}