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

#include "line.hpp"

#include <iostream>

namespace mocc {
real_t Line::distance_to_surface(Point2 p, Direction dir, bool coincident) const
{

    // There can be only one intersection with a line, so if we are already
    // coincident, return infinity
    if (coincident) {
        return std::numeric_limits<real_t>::max();
    }

    // Cast the line into the general form. This might be worth doing ahead of
    // time in the Line struct itself
    real_t a = p1.y - p2.y;
    real_t b = p2.x - p1.x;
    real_t c = p1.x * p2.y - p2.x * p1.y;

    // Evaluate distance to surface
    real_t f = a * p.x + b * p.y + c;

    // Project from location to line
    real_t proj = dir.ox * a + dir.oy * b;

    // Check for point laying on line
    if (std::abs(proj) < REAL_FUZZ) {
        return std::numeric_limits<real_t>::max();
    }

    real_t d = -f / proj;
    if ((d > 0.0) && std::abs(f) > REAL_FUZZ) {
        return d;
    }

    return std::numeric_limits<real_t>::max();

} // Line::distance_to_surface()
} // namespace mocc
