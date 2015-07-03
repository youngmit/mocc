#pragma once


#include "global_config.hpp"

namespace mocc {
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
    };
    
    Angle ToOctant( Angle in, int octant );
}
