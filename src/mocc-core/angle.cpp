#include "angle.hpp"

#include <cassert>
#include <cmath>

#include "global_config.hpp"
#include "constants.hpp"

namespace mocc {
    // Return a new Angle, reflected into the requested octant
    Angle ToOctant( Angle in, int octant ) {
        assert( (0 < octant) & (octant < 9) );

        switch( octant ) {
            case 1:
                return Angle(  fabs(in.ox),
                               fabs(in.oy),
                               fabs(in.oz),
                               in.weight );
            case 2:
                return Angle( -fabs(in.ox),
                               fabs(in.oy),
                               fabs(in.oz),
                               in.weight );
            case 3:
                return Angle( -fabs(in.ox),
                              -fabs(in.oy),
                               fabs(in.oz),
                               in.weight );
            case 4:
                return Angle(  fabs(in.ox),
                              -fabs(in.oy),
                               fabs(in.oz),
                               in.weight );
            case 5:
                return Angle(  fabs(in.ox),
                               fabs(in.oy),
                              -fabs(in.oz),
                               in.weight );
            case 6:
                return Angle( -fabs(in.ox),
                               fabs(in.oy),
                              -fabs(in.oz),
                               in.weight );
            case 7:
                return Angle( -fabs(in.ox),
                              -fabs(in.oy),
                              -fabs(in.oz),
                               in.weight );
            case 8:
                return Angle(  fabs(in.ox),
                              -fabs(in.oy),
                              -fabs(in.oz),
                               in.weight );
        }
        return Angle(0.0, 0.0, 0.0);
    }

    // Change alpha and update ox and oy accordingly
    Angle ModifyAlpha(Angle in, real_t new_alpha) {
        return Angle(new_alpha, in.theta, in.weight);
    }
}
