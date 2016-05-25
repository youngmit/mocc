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

#include <iostream>

#include "pugixml.hpp"

#include "core/pin_mesh.hpp"

using namespace mocc;

TEST(test_cyl) {
    std::string xml_input =
        "<mesh type=\"cyl\" id=\"1\"  pitch=\"1.26\"><radii>0.54 "
        "0.62</radii><sub_radii>4 2</sub_radii><sub_azi>8</sub_azi>";

    pugi::xml_document xml;
    xml.load_string(xml_input.c_str());

    auto pm = PinMeshFactory(xml.child("mesh"));

    std::cout << *pm << std::endl;

    // Check points coincident with interior surfaces.
    // on origin
    CHECK_EQUAL(0, pm->find_reg(Point2(0.0, 0.0), Direction(0.01, HPI)));
    // on positive x-axis
    CHECK_EQUAL(7,
                pm->find_reg(Point2(0.01, 0.0), Direction(TWOPI - 0.01, HPI)));

    // on positive x-axis and outermost ring
    CHECK_EQUAL(
        48, pm->find_reg(Point2(0.62, 0.0), Direction(1.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(
        40, pm->find_reg(Point2(0.62, 0.0), Direction(3.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(
        47, pm->find_reg(Point2(0.62, 0.0), Direction(5.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(
        55, pm->find_reg(Point2(0.62, 0.0), Direction(7.0 * PI / 4.0, HPI)));
}

int main() { return UnitTest::RunAllTests(); }
