#include "pin_mesh_rect.hpp"

#include <iostream>

#include "error.hpp"

namespace mocc {
    PinMesh_Rect::PinMesh_Rect( const pugi::xml_node &input ):
        PinMesh( input ) 
    {
        // Parse the number of x and y divisions
        int ndiv_x = input.child("sub_x").text().as_int(0);
        if (ndiv_x < 1) {
            Error("Failed to read valid number of X divisions in rect pin mesh.");
        }
        int ndiv_y = input.child("sub_y").text().as_int(0);
        if (ndiv_y < 1) {
            Error("Failed to read valid number of Y divisions in rect pin mesh.");
        }

        n_xsreg_ = ndiv_x * ndiv_y;
        n_reg_   = ndiv_x * ndiv_y;

        float dx = pitch_x_/ndiv_x;
        float dy = pitch_y_/ndiv_y;

        for (int i=1; i<ndiv_x; i++) {
            hx_.push_back(i*dx);
        }
        for (int i=1; i<ndiv_y; i++) {
            hy_.push_back(i*dy);
        }
         
    	return;
    }
}
