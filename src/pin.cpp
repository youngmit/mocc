#include "pin.hpp"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include "pugixml.hpp"

#include "error.hpp"
#include "global_config.hpp"
#include "string_utils.hpp"

namespace mocc {
    Pin::Pin( const pugi::xml_node &input, 
            const std::map<int, UP_PinMesh_t> &meshes ):
        id_( input.attribute("id").as_int() )
    { 
        // Check ID validity
        if ( id_ == -1 ) {
            Error( "Failed to read pin ID." );
        }

        // Get pin mesh ID
        mesh_id_ = input.attribute( "mesh" ).as_int( -1 );
        if ( mesh_id_ == -1) {
            Error( "Failed to read pin mesh ID." );
        }
        if ( meshes.count( mesh_id_ ) == 0 ) {
            Error( "Invalid pin mesh ID." );
        }

        pin_mesh_ = meshes.at(mesh_id_).get();

        // Get material IDs
        std::string mats_in = input.child_value();
        std::stringstream inBuf( trim( mats_in ) );
        int mat_id;

        while( !inBuf.eof() ){
            mat_id = -1;
            inBuf >> mat_id;

            mat_IDs_.push_back( mat_id );
        }
        if( inBuf.fail() ){
            Error( "Trouble reading material IDs in pin definition." );
        }

        if ( mat_IDs_.size() != pin_mesh_->n_xsreg() ) {
            Error( "Wrong number of materials specified in pin definition" );
        }
    }
}
