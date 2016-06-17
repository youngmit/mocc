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

#include "circle.hpp"

#include "points.hpp"

#include <limits>

namespace mocc {
real_t Circle::distance_to_surface(Point2 p, Direction dir) const
{
    std::numeric_limits<real_t> lim;

    real_t a = 1.0 - dir.oz * dir.oz;

    if (a == 0.0) {
        return lim.max();
    }

    real_t x = p.x - c.x;
    real_t y = p.y - c.y;

    real_t k   = x * dir.ox + y * dir.oy;
    real_t c   = x * x + y * y - r * r;
    real_t det = k * k - a * c;

    if (det < 0.0) {
        return lim.max();
    }

    // if c ~= 0.0, we are on the surface of the circle. On surfaces, we
    // determine sense w.r.t. direction of travel; the particle is assumed
    // on the side in the direction of travel.
    if (std::abs(c) < REAL_FUZZ) {
        if (k >= 0.0) {
            return lim.max();
        }
        else {
            return (-k + std::sqrt(det)) / a;
        }
    }

    if (c < 0.0) {
        // inside the circle
        return (-k + std::sqrt(det)) / a;
    }
    else {
        // outside the circle
        real_t d = (-k - std::sqrt(det)) / a;
        return d >= 0.0 ? d : lim.max();
    }
    return lim.max();
} // Circle::distance_to_surface()
} // namespace mocc
