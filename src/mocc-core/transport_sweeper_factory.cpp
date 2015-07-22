#include "transport_sweeper_factory.hpp"

namespace mocc {
    TransportSweeper* TransportSweeperFactory( const pugi::xml_node &input,
            const CoreMesh &mesh ) {
        // For now, just make an MoC sweeper. Eventually we will get to PS-style
        // sweepers.

        TransportSweeper* ts = new MoCSweeper( input.child("sweeper"), mesh );

        return ts;
    }
}
