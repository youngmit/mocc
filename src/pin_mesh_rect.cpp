#include "pin_mesh_rect.hpp"

namespace mocc {
    PinMesh_Rect::PinMesh_Rect( const pugi::xml_node &input ):
        PinMesh( input ) {

        n_xsreg_ = 1;
        n_reg_   = 1;
         
    	return;
    }
}
