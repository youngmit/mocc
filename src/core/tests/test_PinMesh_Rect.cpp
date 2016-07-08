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

TEST(test_rect)
{
    std::string xml_input =
        "<mesh type=\"rect\" id=\"1\"  pitch=\"1.26\"><sub_x>5</sub_x><sub_y>"
        "4</sub_y>";

    pugi::xml_document xml;
    xml.load_string(xml_input.c_str());

    auto pm = PinMeshFactory(xml.child("mesh"));

    CHECK(pm);

    std::cout << *pm << std::endl;

    // Check the origin. Should be on a y-normal
    CHECK_EQUAL(12, pm->find_reg(Point2(0.0, 0.0), Direction(PI / 4.0, HPI)));
    CHECK_EQUAL(7,
                pm->find_reg(Point2(0.0, 0.0), Direction(5.0 * PI / 4.0, HPI)));

    // Check somewhere on an X-normal
    CHECK_EQUAL(
        19, pm->find_reg(Point2(0.378, 0.4), Direction(1.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(
        18, pm->find_reg(Point2(0.378, 0.4), Direction(3.0 * PI / 4.0, HPI)));

    // Check somewhere on a corner
    CHECK_EQUAL(6, pm->find_reg(Point2(-0.378, -0.315),
                                Direction(1.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(5, pm->find_reg(Point2(-0.378, -0.315),
                                Direction(3.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(0, pm->find_reg(Point2(-0.378, -0.315),
                                Direction(5.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(1, pm->find_reg(Point2(-0.378, -0.315),
                                Direction(7.0 * PI / 4.0, HPI)));

    // Check along the edges
    // on top edge, pointing out
    CHECK_EQUAL(-1, pm->find_reg(Point2(0.0, 0.63), Direction()));
    // on top edge, pointing in
    CHECK_EQUAL(
        17, pm->find_reg(Point2(0.0, 0.63), Direction(5.0 * PI / 4.0, HPI)));

    // on bottom edge, pointing in
    CHECK_EQUAL(3, pm->find_reg(Point2(0.252, -0.63), Direction()));
    // on bottom edge, pointing out
    CHECK_EQUAL(
        -1, pm->find_reg(Point2(0.252, -0.63), Direction(5.0 * PI / 4.0, HPI)));

    // on right edge, pointing out
    CHECK_EQUAL(-1, pm->find_reg(Point2(0.63, 0.0), Direction()));
    // on bottom edge, pointing in
    CHECK_EQUAL(
        14, pm->find_reg(Point2(0.63, 0.0), Direction(3.0 * PI / 4.0, HPI)));
    CHECK_EQUAL(
        9, pm->find_reg(Point2(0.63, 0.0), Direction(5.0 * PI / 4.0, HPI)));
}

TEST(test_fine_mesh)
{
    std::string xml_input =
        "<mesh id=\"1\" type=\"rect\" pitch=\"10\">"
        "        <sub_x>80</sub_x>"
        "        <sub_y>80</sub_y>"
        "</mesh>";
    pugi::xml_document xml;
    xml.load_string(xml_input.c_str());

    auto pm = PinMeshFactory(xml.child("mesh"));

    CHECK(pm);

    std::cout << *pm << std::endl;

    // two points for legitimate region index assertion test
    CHECK_EQUAL(50,pm->find_reg(Point2(1.250000000000004,-4.9999999999999432)));
    CHECK_EQUAL(6319,pm->find_reg(Point2(4.9375,4.8456249999999992)));
}


int main()
{
    return UnitTest::RunAllTests();
}
