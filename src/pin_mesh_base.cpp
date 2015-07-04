#include "pin_mesh_base.hpp"

#include <iostream>
#include <sstream>

#include "error.hpp"

using std::stringstream;

namespace mocc {
    PinMesh::PinMesh( const pugi::xml_node &input ) {
        // Extract pin id
    	{
    		stringstream inBuf(input.attribute( "id" ).value());
    		inBuf >> id_;
    
    		if(inBuf.fail()) {
    			Error( "Failed to read pin ID." );
    		}
    		if(!inBuf.eof()) {
    			Warn( "Dangling data after pin ID." );
    		}
    	}
    	
    	// Extract pitch
    	{
    		stringstream inBuf( input.attribute( "pitch" ).value() );
    		inBuf >> pitch_x_;
    		// Just treat square pitch for now
    		pitch_y_ = pitch_x_;
    		if( inBuf.fail() ) {
    			Error( "Failed to read pin pitch." );
    		}
    		if( !inBuf.eof() ) {
    			Warn( "Dangling data after pin pitch." );
    		}	
    	}

        return;
    }
}
