#include "angle.hpp"

#include <cassert>
#include <cmath>

#include "global_config.hpp"
#include "constants.hpp"

namespace mocc {
    // Return a new Angle, reflected into the requested octant
    Angle ToOctant( Angle in, int octant ) {
        assert( 0 < octant & octant < 9 );
        Angle out;

        switch( octant ) {
            case 1:
                out.ox = fabs(in.ox);
                out.oy = fabs(in.oy);
                out.oz = fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = acos(out.ox/sin(out.theta));
                break;
            case 2:
                out.ox = -fabs(in.ox);
                out.oy = fabs(in.oy);
                out.oz = fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = acos(out.ox/sin(out.theta));
                break;
            case 3:
                out.ox = -fabs(in.ox);
                out.oy = -fabs(in.oy);
                out.oz = fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = TWOPI-acos(out.ox/sin(out.theta));
                break;
            case 4:
                out.ox = fabs(in.ox);
                out.oy = -fabs(in.oy);
                out.oz = fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = TWOPI-acos(out.ox/sin(out.theta));
                break;
            case 5:
                out.ox = fabs(in.ox);
                out.oy = fabs(in.oy);
                out.oz = -fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = acos(out.ox/sin(out.theta));
                break;
            case 6:
                out.ox = -fabs(in.ox);
                out.oy = fabs(in.oy);
                out.oz = -fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = acos(out.ox/sin(out.theta));
                break;
            case 7:
                out.ox = -fabs(in.ox);
                out.oy = -fabs(in.oy);
                out.oz = -fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = TWOPI-acos(out.ox/sin(out.theta));
                break;
            case 8:
                out.ox = fabs(in.ox);
                out.oy = -fabs(in.oy);
                out.oz = -fabs(in.oz);
                out.theta = acos(out.oz);
                out.alpha = TWOPI-acos(out.ox/sin(out.theta));
                break;
        }
        return out;
    }
}
