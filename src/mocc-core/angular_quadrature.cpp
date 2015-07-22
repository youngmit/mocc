#include "angular_quadrature.hpp"

#include <vector>
#include <string>
#include <iostream>

#include "error.hpp"
#include "level_symmetric.hpp"


namespace mocc {
    AngularQuadrature::AngularQuadrature( const pugi::xml_node &input ) {
        // Make sure we got input
        if( input.empty() ) {
            Error("No input provided for angular quadrature.");
        }

        // Extract the quadrature type
        std::string type_str = input.attribute("type").value();
        if (type_str == "ls") {
            type_ = SN;

            // extract the quadrature order
            int order = input.attribute("order").as_int(-1);

            // Generate angles for octant 1
            angles_ = GenSn( order );
        } else {
            std::cout << input.attribute("type").value() << std::endl;
            Error("Unrecognized angular quadrature type specified.");
        }

        // Store the number of angles per octant
        ndir_oct_ = angles_.size();

        // Expand angles to other octants
        for ( int ioct=2; ioct<=8; ioct++ ) {
            for ( int iang=0; iang < ndir_oct_; iang++ ) {
                Angle a = angles_[iang];
                angles_.push_back( ToOctant(a, ioct) );
            }
        }

        return;
    }

    void AngularQuadrature::modify_angle(int iang, Angle ang ) {
        for ( int ioct=0; ioct<8; ioct++ ){
            angles_[iang + ioct*ndir_oct_] = ToOctant(ang, ioct+1);
        }
    }
    
}
