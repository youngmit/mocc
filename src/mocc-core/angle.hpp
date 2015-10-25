#pragma once

#include <cmath>
#include <iostream>

#include "global_config.hpp"
#include "constants.hpp"

namespace mocc {
    inline real_t RadToDeg( real_t rad ) {
        return 180*(rad*RPI);
    }

    struct Angle {
        /// x-component of the angle
        real_t ox;
        /// y-component of the angle
        real_t oy;
        /// z-component of the angle
        real_t oz;
        /// azimuthal angle
        real_t alpha;
        /// polar cosine
        real_t theta;
        /// quadrature weight
        real_t weight;
        /// Reciprocal of the sine of the polar angle. This is useful for
        /// computing true ray segment length from 2D projected length.
        real_t rsintheta;

        /**
         * Default constructor makes a nonsense angle. Watch out.
         */
        Angle() {}

        /** 
         * Construct using alpha/theta
         */
        Angle( real_t alpha, real_t theta, real_t weight ):
            alpha(alpha),
            theta(theta),
            weight(weight)
        {
            ox = sin(theta)*cos(alpha);
            oy = sin(theta)*sin(alpha);
            oz = cos(theta);
        }

        /**
         * Construct using direction cosines
         */
        Angle( real_t ox, real_t oy, real_t oz, real_t weight):
            ox(ox), oy(oy), oz(oz), weight(weight)
        {
            theta = acos(oz);
            alpha = acos(ox/sin(theta));
            rsintheta = 1.0/sin(theta);
        }
        
        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, const Angle &ang );
    };
    
    Angle ToOctant( Angle in, int octant );

    Angle ModifyAlpha ( Angle in, real_t new_alpha );

    
}
