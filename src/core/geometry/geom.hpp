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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iosfwd>

#include "angle.hpp"
#include "box.hpp"
#include "circle.hpp"
#include "line.hpp"
#include "points.hpp"

#include "core/fp_utils.hpp"
#include "core/global_config.hpp"


namespace {
template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
}

namespace mocc {

// Intersection between a circle and a line segment. Return value carries
// the number of valid intersections that were found
// http://mathworld.wolfram.com/Circle-LineIntersection.html
inline int Intersect(Line l, Circle circ, Point2 &p1, Point2 &p2) {
    int ret = 0;
    real_t u1 = l.p2.x - l.p1.x;
    real_t u2 = l.p2.y - l.p1.y;
    real_t w1 = l.p1.x - circ.c.x;
    real_t w2 = l.p1.y - circ.c.y;

    real_t b = w1 * u1 + w2 * u2;
    real_t c = w1 * w1 + w2 * w2 - circ.r * circ.r;
    if ((c > 0.0) & (b > 0.0)) {
        // no intersection
        p1.ok = false;
        p2.ok = false;
        return 0;
    }

    real_t a = u1 * u1 + u2 * u2;
    real_t discriminant = b * b - a * c;
    if (discriminant < 0.0) {
        // No intersection
        p1.ok = false;
        p2.ok = false;
        return 0;

    } else if (fp_equiv_rel(discriminant, 0.0)) {
        // Tangent. Dont bother
        p1.ok = false;
        p2.ok = false;
        return 0;
    } else {
        // Whoopie, we have a couple of points
        p1.ok = true;
        p2.ok = true;

        real_t ra = 1.0 / a;
        discriminant = sqrt(discriminant);
        real_t t1 = (-b - discriminant) * ra;
        real_t t2 = (-b + discriminant) * ra;
        if ((0.0 < t1) & (t1 < 1.0)) {
            p1 = l.p1;
            p1.x += u1 * t1;
            p1.y += u2 * t1;
            p1.ok = true;
            ret++;
        }
        if ((0.0 < t2) & (t2 < 1.0)) {
            // Make sure that the first point is always valid
            if (ret == 0) {
                assert(false);
                p1 = l.p1;
                p1.x += u1 * t2;
                p1.y += u2 * t2;
                p1.ok = true;
            } else {
                p2 = l.p1;
                p2.x += u1 * t2;
                p2.y += u2 * t2;
                p2.ok = true;
            }
            ret++;
        }
        return ret;
    }
}

inline int Intersect(Line l1, Line l2, Point2 &p) {
    real_t u1 = l1.p2.x - l1.p1.x;
    real_t u2 = l1.p2.y - l1.p1.y;
    real_t v1 = l2.p2.x - l2.p1.x;
    real_t v2 = l2.p2.y - l2.p1.y;
    real_t w1 = l1.p1.x - l2.p1.x;
    real_t w2 = l1.p1.y - l2.p1.y;

    real_t d = u1 * v2 - u2 * v1;
    real_t s = v1 * w2 - v2 * w1;
    real_t t = u1 * w2 - u2 * w1;

    if (fabs(d) < 0.0) {
        // Parallel lines
        p.ok = false;
        return 0;
    }

    s = s / d;
    t = t / d;

    if ((0.0 <= s) & (s <= 1.0) & (0.0 <= t) & (t <= 1.0)) {
        // success
        p = l1.p1;
        p.x += s * u1;
        p.y += s * u2;
        p.ok = true;
        return 1;
    } else {
        // Intersection beyond bounds of line segments
        p.ok = false;
        return 0;
    }
}
}
