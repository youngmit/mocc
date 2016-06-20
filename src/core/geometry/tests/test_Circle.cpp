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

#include "UnitTest++/UnitTest++.h"

#include <cmath>
#include <iostream>
#include <limits>

#include "core/geometry/circle.hpp"

using namespace mocc;

TEST(testCircle)
{
    std::numeric_limits<real_t> lim;
    Circle c(Point2(0.5, 0.5), 0.75);

    // Check coincident point pointing out. Should be max()
    CHECK_EQUAL(lim.max(), c.distance_to_surface(Point2(-0.25, 0.5),
                                                 Direction(-1.0, 0.0, 0.0)));
    CHECK_EQUAL(lim.max(),
                c.distance_to_surface(Point2(-0.25, 0.5),
                                      Direction(-1.0, 0.0, 0.0), true));

    // Check coincident point pointing in and straight across. Should be the
    // diameter
    CHECK_CLOSE(1.5, c.distance_to_surface(Point2(-0.25, 0.5),
                                           Direction(1.0, 0.0, 0.0)),
                REAL_FUZZ);
    // similar to above, but slightly outside
    CHECK_CLOSE(1.5, c.distance_to_surface(Point2(-0.25-0.00000005, 0.5),
                                           Direction(1.0, 0.0, 0.0), false),
                1e-12);

    // Interior point pointing anywhere in the plane
    CHECK_CLOSE(0.75,
                c.distance_to_surface(Point2(0.5, 0.5), Direction(1.0, HPI)),
                REAL_FUZZ);

    // Interior point pointing anywhere out of the plane. Should be the
    // radius/sin(theta)
    CHECK_CLOSE(
        0.75 / std::sin(0.5 * HPI),
        c.distance_to_surface(Point2(0.5, 0.5), Direction(1.0, 0.5 * HPI)),
        REAL_FUZZ);

    // Exterior glancing, should be max()
    CHECK_EQUAL(lim.max(),
                c.distance_to_surface(Point2(-0.25, 0.0),
                                      Direction(HPI + 0.0000001, HPI)));

    // Exterior good
    CHECK_CLOSE(0.15138781886599732327980531686762,
                c.distance_to_surface(Point2(-0.25, 0.0),
                                      Direction(0.5880026035475675, HPI)),
                REAL_FUZZ);

    // Interior, out of the plane. Should be the radius/sin(theta)
    CHECK_CLOSE(
        1.0606601717798212866012665431573,
        c.distance_to_surface(Point2(0.5, 0.5), Direction(0.5, 0.5 * HPI)),
        REAL_FUZZ);
}

int main()
{
    return UnitTest::RunAllTests();
}
