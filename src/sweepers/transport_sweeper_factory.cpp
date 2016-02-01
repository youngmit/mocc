#include "transport_sweeper_factory.hpp"

#include <string>

#include "core/error.hpp"
#include "core/mesh.hpp"

#include "moc/moc_sweeper.hpp"
#include "sn/sn_sweeper_factory.hpp"
#include "cmdo/plane_sweeper_2d3d.hpp"

namespace mocc {
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
            const CoreMesh& mesh ) {
        // Check the input XML for which type of sweeper to make
        std::string type = input.child("sweeper").attribute("type").value();
        if( type == "moc" ) {
            UP_Sweeper_t ts( new moc::MoCSweeper( input.child("sweeper"),
                        mesh ) );
            return ts;
        } else if ( type == "sn" ) {
            auto snts = SnSweeperFactory( input.child("sweeper"), mesh );
            UP_Sweeper_t ts( std::move(snts) );
            return ts;
        } else if ( type == "2d3d" ) {
            UP_Sweeper_t ts( new cmdo::PlaneSweeper_2D3D( input.child("sweeper"),
                        mesh) );
            return ts;
        } else if ( type == "moc_2d3d" ) {
            // Create a 2D3D MoC sweeper by itself. This is really only useful
            // for one-way coupling
            cmdo::MoCSweeper_2D3D *sweeper = new cmdo::MoCSweeper_2D3D(
                    input.child("sweeper"), mesh );
            sweeper->set_self_coupling();
            UP_Sweeper_t ts( sweeper );
            return ts;
        }
        else {
            throw EXCEPT("Failed to detect a valid sweeper type.");
        }
    }
}
