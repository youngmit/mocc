// This exists as the entry point into all pin mesh types, also providing a 
// factory method for generating deferred-type pin mesh objects 

#include "pin_mesh.hpp"

#include <sstream>
#include <iostream>

#include "error.hpp"

namespace mocc {
    // Determine which type of pin to create from an XML object, produce a mesh
    // of the appropriate type and return a shared pointer to the object.
    PinMesh* PinMeshFactory( const pugi::xml_node &input ) {
    	PinMesh* pm;
    	
    	// Extract the type of mesh to make
    	std::string type = input.attribute( "type" ).value();
    
    	if( type == "cyl" ){
    		std::cout << "Found cyl mesh!" << std::endl;
    		pm = new PinMesh_Cyl( input );
    	} else if( type == "rect" ) {
    		std::cout << "Found rect mesh!" << std::endl;
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
}
