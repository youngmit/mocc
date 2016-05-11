/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "pin.hpp"

#include <map>
#include <string>
#include <iostream>

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"
#include "global_config.hpp"
#include "pin_mesh.hpp"
#include "string_utils.hpp"

namespace mocc {
    /**
     * This makes a new Pin from an XML node. By default, a Pin is considered
     * fuel if it has ANY fissile material. The "fuel=\"t|f\"" attribute can be
     * used to override this behavior (useful for non-fuel fissile regions, such
     * as fission chambers and the like).
     */
    Pin::Pin( const pugi::xml_node &input,
            const PinMesh_Map_t &meshes,
            const MaterialLib &mat_lib ) :
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

        // Make sure all of the material IDs are valid, and check for fissile
        // material
        is_fuel_ = false;
        for( auto mid: mat_IDs_ ) {
            if( !mat_lib.has(mid) ) {
                throw EXCEPT("Invalid material specified in pin");
            }
            const auto &mat = mat_lib[mid];
            if( mat.is_fissile() ) {
                is_fuel_ = true;
            }
        }

        // Let user input override the default is_fuel_ setting
        if( !input.attribute("fuel").empty() ) {
            is_fuel_ = input.attribute("fuel").as_bool();
        }
    }

    std::map<int, UP_Pin_t> ParsePins( const pugi::xml_node &input,
            const PinMesh_Map_t &meshes,
            const MaterialLib &mat_lib )
    {
        std::map<int, UP_Pin_t> pins;
        for ( auto pin = input.child( "pin" ); pin;
                pin = pin.next_sibling( "pin" ) )
        {
            // Construct the pin and add to the map
            UP_Pin_t pin_p( new Pin( pin, meshes, mat_lib ) );
            int id = pin_p->id();
            if( pins.find(id) != pins.end() ) {
                std::stringstream msg;
                msg << "Duplicate pin ID (" << id << ")";
                throw EXCEPT(msg.str());
            }
            pins.emplace( id, std::move(pin_p) );
            LogFile << "Pin ID " << id << " done" << std::endl;
        }
        return pins;
    }
}
