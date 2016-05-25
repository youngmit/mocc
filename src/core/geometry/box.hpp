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

#include <array>

#include "core/global_config.hpp"

#include "angle.hpp"
#include "line.hpp"
#include "points.hpp"

namespace mocc {
class Box {
   public:
    /**
     * \note The Lines should be ordered in the same way as the \ref Surface
     * enumeration in \ref constants.hpp
     */
    Box(Point2 p1, Point2 p2)
        : p1_(std::min(p1.x, p2.x), std::min(p1.y, p2.y)),
          p2_(std::max(p1.x, p2.x), std::max(p1.y, p2.y)),
          lines_({Line(Point2(p2_.x, p1_.y), p2_),
                  Line(Point2(p1_.x, p2_.y), p2_),
                  Line(p1_, Point2(p1_.x, p2_.y)),
                  Line(p1_, Point2(p2_.x, p1_.y))}) {
        return;
    }

    Point2 intersect(Point2 p, Angle ang) const;

    std::pair<real_t, Surface> distance_to_surface(Point2 p,
                                                   Direction dir) const;

   private:
    // corners of the box
    Point2 p1_, p2_;
    std::array<Line, 4> lines_;
};
}  // namespace mocc
