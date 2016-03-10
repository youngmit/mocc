#include "pin_mesh_base.hpp"

#include <iostream>

#include "error.hpp"

namespace mocc {
    PinMesh::PinMesh( const pugi::xml_node &input ) {
        // Extract pin id
        id_ = input.attribute( "id" ).as_int(-1);

        if(id_ < 0) {
            throw EXCEPT( "Failed to read pin ID." );
        }

        // Extract pitch
        pitch_x_ = input.attribute( "pitch" ).as_double(-1.0);
        // Just treat square pitch for now
        pitch_y_ = pitch_x_;
        if( pitch_x_ <= 0.0 ) {
            throw EXCEPT( "Failed to read valid pin pitch." );
        }

        return;
    }

    void PinMesh::print( std::ostream &os ) const {
        os << "ID: " << id_ << std::endl;
        os << "X Pitch: " << pitch_x_ << std::endl;
        os << "Y Pitch: " << pitch_y_ << std::endl;
        os << "# of Regions: " << n_reg_ << std::endl;
        os << "# of XS Regions: " << n_xsreg_;
        return;
    }
}
