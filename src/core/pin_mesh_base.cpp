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
