#include "transport_sweeper_factory.hpp"

#include <string>

#include "error.hpp"
#include "mesh.hpp"
#include "moc_sweeper.hpp"
#include "sn_sweeper.hpp"
#include "plane_sweeper_2d3d.hpp"

namespace mocc {
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
            const CoreMesh& mesh ) {
        // Check the input XML for which type of sweeper to make
        std::string type = input.child("sweeper").attribute("type").value();
        if( type == "moc" ) {
            UP_Sweeper_t ts( new MoCSweeper( input.child("sweeper"), mesh ) );
            return ts;
        } else if ( type == "sn" ) {
            UP_Sweeper_t ts( new SnSweeper( input.child("sweeper"), mesh ) );
            return ts;
        } else if ( type == "2d3d" ) {
            UP_Sweeper_t ts( new PlaneSweeper_2D3D( input.child("sweeper"),
                        mesh) ); 
            return ts;
        } else {
            throw EXCEPT("Failed to detect a valid sweeper type.");
        }
    }
}
