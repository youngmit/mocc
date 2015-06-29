#include "transport_sweeper_factory.hpp"

namespace mocc {
    TransportSweeper* TransportSweeperFactory(pugi::xml_node &input){
        // For now, just make an MoC sweeper. Eventually we will get to PS-style
        // sweepers.
        TransportSweeper* ts = new MoCSweeper(input);

        return ts;
    }
}
