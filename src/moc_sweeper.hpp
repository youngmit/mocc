#pragma once

#include "pugixml.hpp"

#include "transport_sweeper.hpp"

namespace mocc {
class MoCSweeper: public TransportSweeper{
public:
    MoCSweeper( const pugi::xml_node &input );
    void sweep(int group);
};
}
