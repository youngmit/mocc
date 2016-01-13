// This exists as the entry point into all pin mesh types, also providing a
// factory method for generating deferred-type pin mesh objects

#include "pin_mesh.hpp"

#include <sstream>
#include <iostream>

#include "error.hpp"
#include "files.hpp"

namespace mocc {
    // Determine which type of pin to create from an XML object, produce a mesh
    // of the appropriate type and return a shared pointer to the object.
    PinMesh* PinMeshFactory( const pugi::xml_node &input ) {
        PinMesh* pm=nullptr;

        // Extract the type of mesh to make
        std::string type = input.attribute( "type" ).value();

        if( type == "cyl" ){
            pm = new PinMesh_Cyl( input );
        } else if( type == "rect" ) {
            pm = new PinMesh_Rect( input );
        } else {
            // I don't recognize the mesh type, error out.
            std::stringstream err;
            err << "Unrecognized mesh type for mesh ID: "
                << input.attribute( "id" ).value();
            Error( err.str().c_str() );
        }

        return pm;
    }

    std::map<int, UP_PinMesh_t> ParsePinMeshes( const pugi::xml_node &input ) {
        std::map<int, UP_PinMesh_t> pin_meshes;
        for (pugi::xml_node mesh = input.child( "mesh" );
                mesh; mesh = mesh.next_sibling( "mesh" ))
        {
            LogFile << "Parsing new pin mesh: ID="
            << mesh.attribute( "id" ).value() << std::endl;
            UP_PinMesh_t pm( PinMeshFactory( mesh ) );
            int id = pm->id();
            pin_meshes.emplace(id, std::move( pm ) );
        }

        return pin_meshes;
    }
}
