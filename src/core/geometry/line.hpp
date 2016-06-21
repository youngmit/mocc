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

#pragma once

#include "util/global_config.hpp"
#include "direction.hpp"
#include "geom_surface.hpp"
#include "points.hpp"

namespace mocc {
struct Line : public GeomSurface {
    Line(Point2 p1, Point2 p2) : p1(p1), p2(p2)
    {
    }

    Point2 p1;
    Point2 p2;

    /**
     * \brief Return the distance to intersection from a point, travelling in a
     * particular direction with the Line
     *
     * \param p a Point2 from which to measure distance.
     * \param dir a Direction along which to measure distance
     * \param coincident whether the given point is to be considered coincident
     * with the line.
     *
     * \todo once things stabilize, document all of the eccentricities
     */
    real_t distance_to_surface(Point2 p, Direction dir,
                               bool coincident=false) const final override;

    friend std::ostream &operator<<(std::ostream &os, Line &l);
};
} // namespace mocc
