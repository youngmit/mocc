#include "angular_quadrature.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "core/error.hpp"
#include "core/level_symmetric.hpp"


namespace mocc {

    const int AngularQuadrature::reflection_[3][8] = { {1,0,3,2,5,4,7,6},
                                                       {3,2,1,0,7,6,5,4},
                                                       {4,5,6,7,0,1,2,3} };

    AngularQuadrature::AngularQuadrature( const pugi::xml_node &input ) {
        // Make sure we got input
        if( input.empty() ) {
            throw EXCEPT("No input provided for angular quadrature.");
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
            std::cout << "'" << type_str << "'" << std::endl;
            throw EXCEPT("Unrecognized angular quadrature type specified.");
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

    std::ostream& operator<<(std::ostream& os,
            const AngularQuadrature &angquad) {
        const int w = 12;
        os << std::setw(w) << "Alpha"
           << std::setw(w) << "Theta"
           << std::setw(w) << "omega x"
           << std::setw(w) << "omega y"
           << std::setw(w) << "omega z" << std::endl;
        for( auto &ang: angquad.angles_ ) {
            os << ang << std::endl;
        }
        return os;
    }
}
