#pragma once

#include "pugixml.hpp"

#include "transport_sweeper.hpp"

namespace mocc {
	/**
	* Peek inside a \<sweeper\> tag to look at the \c type attribute, then
	* generate a TransportSweeeper of the appropriate type using the passed XML
	* node and CoreMesh.
	*/
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
        const CoreMesh& mesh );
}
