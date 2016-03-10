#include "pin.hpp"

#include <map>
#include <string>
#include <iostream>

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"
#include "global_config.hpp"
#include "string_utils.hpp"

namespace mocc {
    Pin::Pin( const pugi::xml_node &input,
            const std::map<int, UP_PinMesh_t> &meshes ):
        id_( input.attribute("id").as_int( 0 ) )
    {
        // Check ID validity
        if ( id_ == 0 ) {
            throw EXCEPT( "Failed to read pin ID." );
        }

        // Get pin mesh ID
        mesh_id_ = input.attribute( "mesh" ).as_int( 0 );
        if ( mesh_id_ == 0 ) {
            throw EXCEPT( "Failed to read pin mesh ID." );
        }
        if ( meshes.count( mesh_id_ ) == 0 ) {
            throw EXCEPT( "Invalid pin mesh ID." );
        }

        pin_mesh_ = meshes.at(mesh_id_).get();

        // Get material IDs
        std::string mats_in = input.child_value();
        mat_IDs_ = explode_string<int>(mats_in);

        if ( (int)mat_IDs_.size() != pin_mesh_->n_xsreg() ) {
            throw EXCEPT( "Wrong number of materials specified in pin "
                    "definition" );
        }
    }

    std::map<int, UP_Pin_t> ParsePins( const pugi::xml_node &input,
            const std::map<int, UP_PinMesh_t> &meshes ) {
        std::map<int, UP_Pin_t> pins;
        for ( auto pin = input.child( "pin" ); pin;
                pin = pin.next_sibling( "pin" ) )
        {
            // Construct the pin and add to the map
            UP_Pin_t pin_p( new Pin( pin, meshes ) );
            int id = pin_p->id();
            pins.emplace( pin_p->id(), std::move(pin_p) );
            LogFile << "Pin ID " << id << " done" << std::endl;
        }
        return pins;
    }
}
