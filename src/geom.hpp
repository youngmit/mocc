#pragma once

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "global_config.hpp"
#include "fp_utils.hpp"
#include "angle.hpp"

#define GEOM_EPS 1e-13

namespace {
    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
}

namespace mocc {
    class Point2 {
    public:
        float_t x;
        float_t y;
        bool ok;

        Point2() {
            ok=false;
            return;
        }
        Point2(float_t x, float_t y):
            x(x), y(y), ok(true)
        {

        }
        
        // return the euclidian distance between the point and another point
        float_t distance(Point2 p) {
            return sqrt((x-p.x)*(x-p.x) + (y-p.y)*(y-p.y));
        }

        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, const Point2 &p) {
            os << "( " << p.x << ", " << p.y << " )";
            return os;
        }

        // Overload the < operator. A point is considered "less than" another
        // point if it has a smaller x-coordinate
        bool operator<(const Point2 &other) const {
            return x < other.x;
        }

        // Point equivalence is based upon approximate floating point
        // arithmetic. Two points are equivalent if their coordinates are "close
        // enough."
        bool operator==(const Point2 &other) const {
            return fp_equiv_abs(x, other.x) & fp_equiv_abs(y, other.y);
        }

        // Vector-like subtraction of two points
        Point2& operator-=(const Point2 &other) {
            x -= other.x;
            y -= other.y;
            return *this;
        }

        // Vector-like subtraction two points
        const Point2 operator-(const Point2 &other) {
            return Point2(*this) -= other;
        }

        // Return the angle, in radians, made by the line from the origin to the
        // point from the positive-x axis. This uses atan2 and some other
        // jiggery to force the result to lie in [0,2PI]
        float_t alpha() const {
            if(y > 0.0) {
                return atan2(y, x);
            } else {
                return atan2(y, x)+TWOPI;
            }
        }
    private:
    };

    class Box {
    public:
        Box( Point2 p1, Point2 p2 ) {
            p1_ = Point2( std::min(p1.x, p2.x), std::min(p1.y, p2.y) );
            p2_ = Point2( std::max(p1.x, p2.x), std::max(p1.y, p2.y) );

            return;
        }

        Point2 intersect( Point2 p, Angle ang ) {
            // Project ox/oy to 2D plane
            float_t ox = cos(ang.alpha);
            float_t oy = sin(ang.alpha);

            // Dont use this code for astronomy stuff
            float_t d_min = 1.0e12;
            float_t d;
            float_t x, y;

            Point2 p_out;

            // Check distance to x-normal planes
            x = p1_.x;
            d = ( p1_.x - p.x ) / ox;
            y = p.y+oy*d;
            if ( (d > GEOM_EPS) & (d < d_min) & (y > p1_.y) & (y < p2_.y )) {
                d_min = d;
                p_out.x = x;
                p_out.y = y;
                p_out.ok = true;
            }

            x = p2_.x;
            d = ( p2_.x - p.x ) / ox;
            y = p.y+oy*d;
            if ( (d > GEOM_EPS) & (d < d_min) & (y > p1_.y) & (y < p2_.y) ) {
                d_min = d;
                p_out.x = x;
                p_out.y = y;
                p_out.ok = true;
            }

            // Check distance to y-normal planes
            y = p1_.y;
            d = ( p1_.y - p.y ) / oy;
            x = p.x+ox*d;
            if ( (d > GEOM_EPS) & (d < d_min) & (x > p1_.x) & (x < p2_.x) ) {
                d_min = d;
                p_out.x = x;
                p_out.y = y;
                p_out.ok = true;
            }

            y = p2_.y;
            d = ( p2_.y - p.y ) / oy;
            x = p.x+ox*d;
            if ( (d > GEOM_EPS) & (d < d_min) & (x > p1_.x) & (x < p2_.x) ) {
                d_min = d;
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

    struct Circle { 
        Point2 c;
        float_t r;
        Circle(Point2 c, float_t r): c(c), r(r) { }
    };

    struct Line {
        Point2 p1;
        Point2 p2;
        Line(Point2 p1, Point2 p2): p1(p1), p2(p2) { }
    };

    inline Point2 Midpoint( const Point2 p1, const Point2 p2 ) {
        return Point2(0.5*(p1.x + p2.x), 0.5*(p1.y + p2.y));
    }

    // Intersection between a circle and a line segment. Return value carries
    // the number of valid intersections that were found 
    // http://mathworld.wolfram.com/Circle-LineIntersection.html 
    inline int Intersect( Line l, Circle circ, Point2 &p1, Point2 &p2 ) {
        int ret = 0;
        float_t u1 = l.p2.x - l.p1.x;
        float_t u2 = l.p2.y - l.p1.y;
        float_t w1 = l.p1.x - circ.c.x;
        float_t w2 = l.p1.y - circ.c.y;

        float_t b = w1*u1 + w2*u2;
        float_t c = w1*w1 + w2*w2 - circ.r*circ.r;
        if ( (c > 0.0) & (b > 0.0) ) {
            // no intersection
            p1.ok = false;
            p2.ok = false;
            return 0;
        }

        float_t a = u1*u1 + u2*u2;
        float_t discriminant = b*b-a*c;
        if ( discriminant < 0.0 ) {
            // No intersection
            p1.ok = false;
            p2.ok = false;
            return 0;
            
        } else if ( fp_equiv_rel(discriminant, 0.0) ) {
            // Tangent. Dont bother
            p1.ok = false;
            p2.ok = false;
            return 0;
        } else {
            // Whoopie, we have a couple of points
            p1.ok = true;
            p2.ok = true;
            
            float_t ra = 1.0/a;
            discriminant = sqrt(discriminant);
            float_t t1 = (-b-discriminant)*ra;
            float_t t2 = (-b+discriminant)*ra;
            if( (0.0 < t1) & (t1 < 1.0) ) {
                p1 = l.p1;
                p1.x += u1*t1;
                p1.y += u2*t1;
                p1.ok = true;
                ret++;
            }
            if( (0.0 < t2) & (t2 < 1.0) ) {
                // Make sure that the first point is always valid
                if (ret == 0) {
                    assert(false);
                    p1 = l.p1;
                    p1.x += u1*t2;
                    p1.y += u2*t2;
                    p1.ok = true;
                } else {
                    p2 = l.p1;
                    p2.x += u1*t2;
                    p2.y += u2*t2;
                    p2.ok = true;
                }
                ret++;
            }
            return ret;
        }
    }


    inline int Intersect( Line l1, Line l2, Point2 &p ) {
        float_t u1 = l1.p2.x - l1.p1.x;
        float_t u2 = l1.p2.y - l1.p1.y;
        float_t v1 = l2.p2.x - l2.p1.x;
        float_t v2 = l2.p2.y - l2.p1.y;
        float_t w1 = l1.p1.x - l2.p1.x;
        float_t w2 = l1.p1.y - l2.p1.y;

        float_t d = u1*v2 - u2*v1;
        float_t s = v1*w2 - v2*w1;
        float_t t = u1*w2 - u2*w1;

        if ( fabs(d) < 0.0 ) {
            // Parallel lines
            p.ok = false;
            return 0;
        }

        s = s/d;
        t = t/d;

        if ( (0.0 <= s) & (s <= 1.0) & (0.0 <= t) & (t <= 1.0) ) {
            // success
            p = l1.p1;
            p.x += s*u1;
            p.y += s*u2;
            p.ok = true;
            return 1;
        } else {
            // Intersection beyond bounds of line segments
            p.ok = false;
            return 0;
        }
    }
}