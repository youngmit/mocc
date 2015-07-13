#pragma once

#include <cmath>
#include <iostream>

#include "global_config.hpp"
#include "constants.hpp"

namespace mocc {
    inline float_t RadToDeg( float_t rad ) {
        return 180*(rad*RPI);
    }

    struct Angle {
        // x- y- z- components of the angle
        float_t ox;
        float_t oy;
        float_t oz;
        // azimuthal angle
        float_t alpha;
        // polar cosine
        float_t theta;
        // quadrature weight
        float_t weight;

        // Construct using alpha/theta
        Angle( float_t alpha, float_t theta, float_t weight ):
            alpha(alpha),
            theta(theta),
            weight(weight)
        {
            ox = sin(theta)*cos(alpha);
            oy = sin(theta)*sin(alpha);
            oz = cos(theta);
        }

        // Construct using direction cosines
        Angle( float_t ox, float_t oy, float_t oz, float_t weight):
            ox(ox), oy(oy), oz(oz), weight(weight)
        {
            theta = acos(oz);
            alpha = acos(ox/sin(theta));
        }
        
        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, const Angle &ang ) {
            os << RadToDeg(ang.alpha) << "\t" << RadToDeg(ang.theta) << "\t"
               << ang.ox << "   \t" << ang.oy << "   \t" << ang.oz << "   \t"
               << ang.weight;
            return os;
        }

    };
    
    Angle ToOctant( Angle in, int octant );

    Angle ModifyAlpha ( Angle in, float_t new_alpha );

    
}
