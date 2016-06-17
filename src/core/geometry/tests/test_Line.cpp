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

#include <limits>

#include "core/geometry/line.hpp"

using namespace mocc;

TEST(test_line)
{
    std::numeric_limits<real_t> lim;
    {
        Line l(Point2(0.0, -10.0), Point2(0.0, 10.0));

        // Simplest case. Hit head on from unit distance
        CHECK_CLOSE(1.0, l.distance_to_surface(Point2(-1.0, 0.0),
                                               Direction(1.0, 0.0, 0.0)),
                    REAL_FUZZ);

        // Reverse direction, should return max float value
        CHECK_EQUAL(lim.max(),
                    l.distance_to_surface(Point2(-1.0, 0.0),
                                          Direction(-1.0, 0.0, 0.0)));

        // A little more complex, hit at 45 degrees, should be root 2
        CHECK_CLOSE(1.4142135623730950488016887242097,
                    l.distance_to_surface(
                        Point2(-1.0, 0.0),
                        Direction(0.70710678118654752440084436210485,
                                  0.70710678118654752440084436210485, 0.0)),
                    REAL_FUZZ);

        // Pretty much the same, Direction defined diferently
        CHECK_CLOSE(1.4142135623730950488016887242097,
                    l.distance_to_surface(
                        Point2(-1.0, 0.0),
                        Direction(0.78539816339744830961566084581988, HPI)),
                    REAL_FUZZ);

        // Add some vertical component
        CHECK_CLOSE(2.0, l.distance_to_surface(
                             Point2(-1.0, 0.0),
                             Direction(0.78539816339744830961566084581988,
                                       0.78539816339744830961566084581988)),
                    REAL_FUZZ);

        // Conincidence
        CHECK_EQUAL(
            lim.max(),
            l.distance_to_surface(Point2(1.0, 5.0), Direction(HPI, 0.5 * HPI)));
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
