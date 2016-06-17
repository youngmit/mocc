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

#include "util/fp_utils.hpp"
#include "util/global_config.hpp"
#include "core/constants.hpp"

namespace mocc {
struct Point2 {
public:
    real_t x;
    real_t y;
    bool ok;

    Point2()
    {
        ok = false;
        return;
    }
    Point2(real_t x, real_t y) : x(x), y(y), ok(true)
    {
    }

    // return the euclidian distance between the point and another point
    real_t distance(Point2 p)
    {
        return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
    }

    // Provide stream insertion support
    friend std::ostream &operator<<(std::ostream &os, const Point2 &p);

    /**
     * A point is considered "less than" another point if it has a smaller
     * y-coordinate.
     */
    bool operator<(const Point2 &other) const
    {
        return y < other.y;
    }

    // Point equivalence is based upon approximate floating point
    // arithmetic. Two points are equivalent if their coordinates are "close
    // enough."
    bool operator==(const Point2 &other) const
    {
        return fp_equiv_abs(x, other.x) && fp_equiv_abs(y, other.y);
    }

    // Vector-like subtraction of two points
    Point2 &operator-=(const Point2 &other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // Vector-like subtraction two points
    const Point2 operator-(const Point2 &other)
    {
        return Point2(*this) -= other;
    }

    // Return the angle, in radians, made by the line from the origin to the
    // point from the positive-x axis. This uses atan2 and some other
    // jiggery to force the result to lie in [0,2PI]
    real_t alpha() const
    {
        if (y > 0.0) {
            return atan2(y, x);
        }
        else {
            return atan2(y, x) + TWOPI;
        }
    }

private:
};

struct Point3 {
public:
    real_t x;
    real_t y;
    real_t z;
    bool ok;

    Point3() : x(0.0), y(0.0), z(0.0), ok(false)
    {
        return;
    }

    Point3(real_t x, real_t y, real_t z) : x(x), y(y), z(z), ok(true)
    {
        return;
    }

    // return the euclidian distance between the point and another point
    real_t distance(Point3 p)
    {
        return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) +
                    (z - p.z) * (z - p.z));
    }

    /**
     * A \ref Point3 is considered "less than" another if it is closer to
     * the origin.
     */
    bool operator<(const Point3 &other) const
    {
        return (x * x + y * y + z * z) <
               (other.x * other.x + other.y * other.y + other.z * other.z);
    }

    // Point equivalence is based upon approximate floating point
    // arithmetic. Two points are equivalent if their coordinates are "close
    // enough."
    bool operator==(const Point3 &other) const
    {
        return fp_equiv_abs(x, other.x) && fp_equiv_abs(y, other.y) &&
               fp_equiv_abs(z, other.z);
    }

    // Vector-like subtraction of two points
    Point3 &operator-=(const Point3 &other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    // Vector-like subtraction two points
    const Point3 operator-(const Point3 &other)
    {
        return Point3(*this) -= other;
    }

    /**
     * \brief Return a \ref Point2 instance containing the 2-dimensional
     * component of the \ref Point3
     */
    Point2 to_2d() const
    {
        return Point2(x, y);
    }
};

inline Point2 Midpoint(const Point2 p1, const Point2 p2)
{
    return Point2(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y));
}

} // namespace mocc
