#pragma once

#include "pugixml.hpp"

#include "mocc-core/transport_sweeper.hpp"

namespace mocc {
    /**
    * Peek inside a \<sweeper\> tag to look at the \c type attribute, then
    * generate a \ref TransportSweeper of the appropriate type using the passed
    * XML node and \ref CoreMesh.
    */
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
        const CoreMesh& mesh );
}
