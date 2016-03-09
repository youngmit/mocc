#include "transport_sweeper_factory.hpp"

#include <string>

#include "core/error.hpp"
#include "core/files.hpp"
#include "core/mesh.hpp"

#include "moc/moc_sweeper.hpp"
#include "sn_sweeper_factory.hpp"
#include "cmdo/plane_sweeper_2d3d.hpp"

namespace mocc {
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
            const CoreMesh& mesh ) {

        LogFile << "Generating transport sweeper..." << std::endl; 

        UP_Sweeper_t sweeper;

        // Check the input XML for which type of sweeper to make
        std::string type = input.child("sweeper").attribute("type").value();
        if( type == "moc" ) {
            LogScreen << "Using an MoC sweeper" << std::endl;
            UP_Sweeper_t ts( new moc::MoCSweeper( input.child("sweeper"),
                        mesh ) );
            sweeper = std::move( ts );
        } else if ( type == "sn" ) {
            LogScreen << "Using an Sn sweeper" << std::endl;
            auto snts = SnSweeperFactory( input.child("sweeper"), mesh );
            sweeper = std::move(snts);
        } else if ( type == "2d3d" ) {
            LogScreen << "Using a 2D3D wweeper" << std::endl;
            UP_Sweeper_t ts( new cmdo::PlaneSweeper_2D3D( input.child("sweeper"),
                        mesh) );
            sweeper = std::move( ts );
        } else if ( type == "moc_2d3d" ) {
            LogScreen << "Using a standalone 2D3D MoC sweeper" << std::endl;
            // Create a 2D3D MoC sweeper by itself. This is really only useful
            // for one-way coupling
            cmdo::MoCSweeper_2D3D *swp = new cmdo::MoCSweeper_2D3D(
                    input.child("sweeper"), mesh );
            swp->set_self_coupling();
            sweeper.reset( swp );
        } else {
            throw EXCEPT("Failed to detect a valid sweeper type.");
        }

        LogFile << "Done generating transport sweeper." << std::endl;

        return sweeper;
    }
}
