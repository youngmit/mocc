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
#include <iostream>
#include "pugixml.hpp"
#include "constants.hpp"
#include "util/error.hpp"
#include "util/global_config.hpp"

namespace mocc {
Angle::Angle(const pugi::xml_node &input)
{
    // All angles need to specify a weight and a direction. The direction
    // may be specified using direction cosines (x, y and z components) or
    // as a polar and an azimuthal angle. If a mixture of the two is use, we
    // will throw an error.

    real_t weight = input.attribute("weight").as_double(0.0);
    if (weight <= 0.0) {
        throw EXCEPT("Invalid angle weight specified.");
    }

    if (!input.attribute("ox").empty()) {
        // Make sure that there arent any polar or azimuthal angles
        if (!input.attribute("alpha").empty() ||
            !input.attribute("theta").empty()) {
            throw EXCEPT("An angle appears to be over-defined (both "
                         "direction cosines and polar/azimuthal angles are "
                         "specified)");
        }
        real_t ox = input.attribute("ox").as_double(2.0);
        real_t oy = input.attribute("oy").as_double(2.0);
        real_t oz = input.attribute("oz").as_double(2.0);
        // Check the direction cosines for validity. (-1,+1) && sqrt(sum of
        // squares) == 1
        // The valid range for direction cosines is exclusive, since we
        // assume that angles lie unambiguously within an octant, and
        // therefore should not lie directly on any of the axes.
        if (ox <= -1.0 || ox >= 1.0) {
            throw EXCEPT("Invalid ox in angle.");
        }
        if (oy <= -1.0 || oy >= 1.0) {
            std::cout << oy << std::endl;
            throw EXCEPT("Invalid oy in angle.");
        }
        if (oz <= -1.0 || oz >= 1.0) {
            throw EXCEPT("Invalid oz in angle.");
        }
        if (!fp_equiv_ulp(sqrt(ox * ox + oy * oy + oz * oz), 1.0)) {
            throw EXCEPT("Direction cosines dont lie on the unit sphere.");
        }

        *this = Angle(ox, oy, oz, weight);
        return;
    }
    else if (!input.attribute("alpha").empty()) {
        // Make sure that there arent any direction cosines
        if (!input.attribute("ox").empty() || !input.attribute("oy").empty() ||
            !input.attribute("oz").empty()) {
            throw EXCEPT("An angle appears to be over-defined (both "
                         "direction cosines and polar/azimuthal angles are "
                         "specified)");
        }

        real_t theta = input.attribute("theta").as_double();
        real_t alpha = input.attribute("alpha").as_double();

        // Check the angles for validity. Theta in (-pi/2, pi/2), alpha in
        // [0, 2*pi)
        // Theta is inclusive, since most 2d sweepers would explode with
        // vertical polar angles. Azimuthal angle is allowed to be zero, but
        // less than 2*pi, to avoid overlapping angles on the positive
        // x-axis.
        if (theta <= -HPI || theta >= HPI) {
            throw EXCEPT("Invalid polar angle.");
        }
        if ((alpha < 0.0) || (alpha >= TWOPI) || (alpha == 0.0) ||
            (alpha == HPI) || (alpha == PI) || (alpha == 3.0 * HPI)) {
            throw EXCEPT("Invalid azimuthal angle.");
        }

        *this = Angle(alpha, theta, weight);
        return;
    }
    else {
        throw EXCEPT("No valid direction specified for angle.");
    }
    return;
}

// Return a new Angle, reflected into the requested octant
Angle Angle::to_octant(int octant) const
{
    return Angle(Direction::to_octant(octant), weight);
}

std::ostream &operator<<(std::ostream &os, const Angle &ang)
{
    const int w = 12;
    os << std::setw(w) << RadToDeg(ang.alpha) << std::setw(w)
       << RadToDeg(ang.theta) << std::setw(w) << ang.ox << std::setw(w)
       << ang.oy << std::setw(w) << ang.oz << std::setw(w) << ang.weight
       << std::setw(w) << ang.rsintheta;
    return os;
}
}
