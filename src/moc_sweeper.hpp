#pragma once

#include "pugixml.hpp"

#include "transport_sweeper.hpp"

namespace mocc {
class MoCSweeper: public TransportSweeper{
public:
    MoCSweeper(pugi::xml_node &input);
    void sweep(int group);
};
}
