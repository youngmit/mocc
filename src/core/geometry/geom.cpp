/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distrealributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "geom.hpp"

#include <iostream>
#include <limits>

namespace mocc {
std::ostream &operator<<(std::ostream &os, const Point2 &p) {
    os << "[ " << p.x << ", " << p.y << " ]";
    return os;
}

std::ostream &operator<<(std::ostream &os, Line &l) {
    os << "[" << l.p1 << ", " << l.p2 << "]";
    return os;
}

real_t Circle::distance_to_surface(Particle p) const {
    real_t d;
    std::numeric_limits<real_t> lim;

    real_t a = 1.0 - p.direction.oz * p.direction.oz;

    if (a == 0.0) {
        return lim.max();
    }

    real_t x = p.location.x - c.x;
    real_t y = p.location.y - c.y;

    real_t k = x * p.direction.ox + y * p.direction.oy;
    real_t c = x * x + y * y - r * r;
    real_t det = k * k - a * c;

    if (det < 0.0) {
        return lim.max();
    }

    // if c ~= 0.0, we are on the surface of the circle. On surfaces, we
    // determine sense w.r.t. direction of travel; the particle is assumed
    // on the side in the direction of travel.
    if (std::abs(c) < 4.0 * lim.epsilon()) {
        if (k >= 0.0) {
            return lim.max();
        } else {
            return (-k + std::sqrt(det)) / a;
        }
    }

    if (c < 0.0) {
        // inside the circle
        return (-k + std::sqrt(det)) / a;
    } else {
        // outside the circle
        real_t d = (-k - std::sqrt(det)) / a;
        return d >= 0.0 ? d : lim.max();
    }
    return lim.max();
} // Circle::distance_to_surface()

Line::distance_to_surface( Particle p ) const {
    const std::numeric_limits<real_t> lim;

    // Cast the line into the general form. This might be worth doing ahead of
    // time in the Line struct itself
    real_t a = p1.y - p2.y;
    real_t b = p2.x - p1.x;
    real_t c = p1.x*p2.y - p2.x*p1.y;

    // Evaluate distance to surface
    real_t f = a*p.location.x + b*p.location.y + c;

    // Project from location to line
    real_t p = p.direction.ox*a + p.direction.oy*b;

    // Check for point laying on line
    if(std::abs(p) < 4.0*lim.epsilon()) {
        return lim.max();
    }

    real_t d = -f/p;
    if( d >= 0.0 ) {
        return d;
    }

    return lim.max();

} // Line::distance_to_surface()
}
