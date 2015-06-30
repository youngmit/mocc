#include "core_mesh.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include "files.hpp"
#include "error.hpp"
#include "string_utils.hpp"
#include "assembly.hpp"


using std::cout;
using std::endl;
using std::stringstream;


namespace mocc {
    CoreMesh::CoreMesh(pugi::xml_node &input) {
        // Parse meshes
        for (pugi::xml_node mesh = input.child( "mesh" );
             mesh; 
             mesh = mesh.next_sibling( "mesh" )) {
            LogFile << "Parsing new pin mesh: ID=" 
            << mesh.attribute( "id" ).value() << endl;
            UP_PinMesh_t pm( PinMeshFactory( mesh ) );
            cout << "pin mesh id: " << pm->id() << endl;
            pin_meshes_.emplace(pm->id(), std::move( pm ) );
        }
        
        // Parse Material Library
        std::string matLibName = 
            input.child( "material_lib" ).attribute( "path" ).value();
        cout << "Found material library specification: " << matLibName << endl;
        FileScrubber matLibFile( matLibName.c_str(), "!" );
        mat_lib_ = MaterialLib( matLibFile );
        
        // Parse material IDs
        for (pugi::xml_node mat = input.child( "material_lib" ).child( "material" ); 
                mat; mat = mat.next_sibling( "material" )){
            cout << mat.attribute( "id" ).value() << " "
                 << mat.attribute( "name" ).value() << endl;
            mat_lib_.assignID( mat.attribute( "id" ).as_int(),
                               mat.attribute( "name" ).value() );
        }
        
        // Parse pins
        for ( pugi::xml_node pin = input.child( "pin" ); pin; 
                pin = pin.next_sibling( "pin" ) ) {
            // Get pin ID
            int pin_id = pin.attribute( "id" ).as_int( -1 );
            if ( pin_id == -1 ) {
                Error( "Failed to read pin ID." );
            }

            // Get pin mesh ID
            int mesh_id = pin.attribute( "mesh" ).as_int( -1 );
            if ( mesh_id == -1) {
                Error( "Failed to read pin mesh ID." );
            }
            if ( pin_meshes_.count( mesh_id ) == 0 ) {
                Error( "Invalid pin mesh ID." );
            }

            // Get material IDs
            std::string mats_in = pin.child_value();
            stringstream inBuf( trim( mats_in ) );
            int mat_id;
            VecI mats;

            while( !inBuf.eof() ){
                mat_id = -1;
                inBuf >> mat_id;

                mats.push_back( mat_id );
            }
            if( inBuf.fail() ){
                Error( "Trouble reading material IDs in pin definition." );
            }
            if ( mats.size() != pin_meshes_[pin_id]->n_xsreg() ) {
                Error( "Wrong number of materials specified in pin definition" );
            }

            // Construct the pin and add to the map
            UP_Pin_t pin_p( new Pin( pin_id, pin_meshes_[mesh_id].get(), mats   ) );
            pins_.emplace( pin_id, std::move(pin_p) );
        }

        // Parse lattices
        for ( pugi::xml_node lat = input.child( "lattice" ); lat;
                lat = input.next_sibling( "lattice" )) {
            Lattice lattice( lat, pins_ );

            lattices_.emplace( lattice.id(), lattice );
        }

        // Parse assemblies
        for ( pugi::xml_node asy = input.child("assembly"); asy;
                asy = input.next_sibling("assembly") ) {
            Assembly assembly( asy, lattices_ );
        }

        // Parse core

        return;
    }

    CoreMesh::~CoreMesh() {
        std::cout << "Destroying CoreMesh" << std::endl;
        return;
    }
}
