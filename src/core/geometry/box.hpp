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

#include "core/global_config.hpp"

#include "points.hpp"

namespace mocc {
class Box {
   public:
    Box(Point2 p1, Point2 p2) {
        p1_ = Point2(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
        p2_ = Point2(std::max(p1.x, p2.x), std::max(p1.y, p2.y));

        return;
    }

    Point2 intersect(Point2 p, Angle ang) {
        // Project ox/oy to 2D plane
        real_t ox = cos(ang.alpha);
        real_t oy = sin(ang.alpha);

        // Dont use this code for astronomy stuff
        real_t d_min = 1.0e12;
        real_t d;
        real_t x, y;

        Point2 p_out;

        // Check distance to x-normal planes
        x = p1_.x;
        d = (p1_.x - p.x) / ox;
        y = p.y + oy * d;
        if ((fabs(d) > GEOM_EPS) && (d < d_min) && (y > p1_.y) && (y < p2_.y)) {
            d_min = fabs(d);
            p_out.x = x;
            p_out.y = y;
            p_out.ok = true;
        }

        x = p2_.x;
        d = (p2_.x - p.x) / ox;
        y = p.y + oy * d;
        if ((fabs(d) > GEOM_EPS) && (d < d_min) && (y > p1_.y) && (y < p2_.y)) {
            d_min = fabs(d);
            p_out.x = x;
            p_out.y = y;
            p_out.ok = true;
        }

        // Check distance to y-normal planes
        y = p1_.y;
        d = (p1_.y - p.y) / oy;
        x = p.x + ox * d;
        if ((fabs(d) > GEOM_EPS) && (d < d_min) && (x > p1_.x) && (x < p2_.x)) {
            d_min = fabs(d);
            p_out.x = x;
            p_out.y = y;
            p_out.ok = true;
        }

        y = p2_.y;
        d = (p2_.y - p.y) / oy;
        x = p.x + ox * d;
        if ((fabs(d) > GEOM_EPS) && (d < d_min) && (x > p1_.x) && (x < p2_.x)) {
            d_min = fabs(d);
            p_out.x = x;
            p_out.y = y;
            p_out.ok = true;
        }

        return p_out;
    }

   private:
    // corners of the box
    Point2 p1_, p2_;
};
} // namespace mocc
