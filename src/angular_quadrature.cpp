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
        }

        // Store the number of angle per octant
        ndir_oct_ = angles_.size();

        // Expand angles to other octants
        for ( int ioct = 2; ioct<=8; ioct++ ) {
            for ( int iang=0; iang < ndir_oct_; iang++ ) {
                Angle a = angles_[iang];
                angles_.push_back( ToOctant(a, ioct) );
            }
        }

        return;
    }
}
