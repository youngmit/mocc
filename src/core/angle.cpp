/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "angle.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>

#include "global_config.hpp"
#include "constants.hpp"

namespace mocc {
    // Return a new Angle, reflected into the requested octant
    Angle Angle::to_octant( int octant ) const {
        assert( (0 < octant) & (octant < 9) );

        switch( octant ) {
            case 1:
                return Angle(  fabs(ox),
                               fabs(oy),
                               fabs(oz),
                               weight );
            case 2:
                return Angle( -fabs(ox),
                               fabs(oy),
                               fabs(oz),
                               weight );
            case 3:
                return Angle( -fabs(ox),
                              -fabs(oy),
                               fabs(oz),
                               weight );
            case 4:
                return Angle(  fabs(ox),
                              -fabs(oy),
                               fabs(oz),
                               weight );
            case 5:
                return Angle(  fabs(ox),
                               fabs(oy),
                              -fabs(oz),
                               weight );
            case 6:
                return Angle( -fabs(ox),
                               fabs(oy),
                              -fabs(oz),
                               weight );
            case 7:
                return Angle( -fabs(ox),
                              -fabs(oy),
                              -fabs(oz),
                               weight );
            case 8:
                return Angle(  fabs(ox),
                              -fabs(oy),
                              -fabs(oz),
                               weight );
        }
        return Angle(0.0, 0.0, 0.0);
    }

    std::ostream& operator<<(std::ostream& os, const Angle &ang ) {
        const int w = 12;
            os << std::setw(w) << RadToDeg(ang.alpha)
               << std::setw(w) << RadToDeg(ang.theta)
               << std::setw(w) << ang.ox
               << std::setw(w) << ang.oy
               << std::setw(w) << ang.oz
               << std::setw(w) << ang.weight
               << std::setw(w) << ang.rsintheta;
            return os;
        }
}
