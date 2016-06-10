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

#include "angle.hpp"
#include "constants.hpp"
#include "error.hpp"

using mocc::Angle;
using mocc::Exception;

TEST(testAngle)
{
    // Make angles with alpha and theta
    Angle a(PI / 6.0, PI / 4.0, 0.564);
    CHECK_CLOSE(0.523598775598299, a.alpha, 0.000000000001);
    CHECK_CLOSE(0.785398163397448, a.theta, 0.000000000001);
    CHECK_CLOSE(0.612372435695794, a.ox, 0.00000000001);
    CHECK_CLOSE(0.353553390593274, a.oy, 0.00000000001);
    CHECK_CLOSE(0.707106781186548, a.oz, 0.00000000001);
    CHECK_EQUAL(0.564, a.weight);

    // Make the same angle with the equivalent direction cosines
    Angle a2(0.612372435695794, 0.353553390593274, 0.707106781186548, 0.564);
    CHECK(a == a2);

    // Reflect into another octant
    Angle a3 = a2.to_octant(4);
    CHECK_CLOSE(5.759586531581287, a3.alpha, 0.000000000001);
    CHECK_CLOSE(0.785398163397448, a3.theta, 0.000000000001);
    CHECK_CLOSE(0.612372435695794, a3.ox, 0.00000000001);
    CHECK_CLOSE(-0.353553390593274, a3.oy, 0.00000000001);
    CHECK_CLOSE(0.707106781186548, a3.oz, 0.00000000001);
    CHECK_EQUAL(0.564, a.weight);
}

TEST(testXMLAngle)
{
    {
        // this should fail
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.0\"/>");
        try {
            Angle a(xml.child("angle"));
            CHECK(false);
        }
        catch (Exception e) {
            // std::cout << e.what() << std::endl;
        }
        catch (...) {
            CHECK(false);
        }
    }
    {
        // this should fail (overdefined)
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.1\" ox=\"0.612372435695794\" "
                        "oy=\"0.353553390593274\" "
                        "oz=\"0.707106781186548\" theta=\"0.5\" />");
        try {
            Angle a(xml.child("angle"));
            throw 0;
        }
        catch (Exception e) {
            std::cout << e.what() << std::endl;
        }
        catch (int) {
            CHECK(false);
        }
    }

    {
        // This should fail. Invalid angle, not on unit sphere
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.1\" ox=\"0.612372435695794\" "
                        "oy=\"0.353553390593274\" "
                        "oz=\"0.708106781186548\" />");

        try {
            Angle a(xml.child("angle"));
            throw 0;
        }
        catch (Exception e) {
            std::cout << e.what() << std::endl;
        }
        catch (int) {
            CHECK(false);
        }
    }
    {
        // This should fail. Forbidden angle, using cosines
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.1\" ox=\"1.0\" "
                        "oy=\"0.0\" "
                        "oz=\"0.0\" />");

        try {
            Angle a(xml.child("angle"));
            throw 0;
        }
        catch (Exception e) {
            std::cout << e.what() << std::endl;
        }
        catch (int) {
            CHECK(false);
        }
    }
    {
        // This should fail. Forbidden angle, using alpha/theta
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.1\" alpha=\"3.1415926535897932\" "
                        "theta=\"0.5\" />");

        try {
            Angle a(xml.child("angle"));
            throw 0;
        }
        catch (Exception e) {
            std::cout << e.what() << std::endl;
        }
        catch (int) {
            CHECK(false);
        }
    }
    {
        // This should pass.
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.1\" ox=\"0.612372435695794\" "
                        "oy=\"0.353553390593274\" "
                        "oz=\"0.707106781186548\" />");
        Angle a(xml.child("angle"));

        CHECK_CLOSE(0.523598775598299, a.alpha, 0.000000000001);
        CHECK_CLOSE(0.785398163397448, a.theta, 0.000000000001);
        CHECK_CLOSE(0.612372435695794, a.ox, 0.00000000001);
        CHECK_CLOSE(0.353553390593274, a.oy, 0.00000000001);
        CHECK_CLOSE(0.707106781186548, a.oz, 0.00000000001);
        CHECK_EQUAL(0.1, a.weight);
    }
    {
        // This should pass.
        pugi::xml_document xml;
        xml.load_string("<angle weight=\"0.1\" alpha=\"0.523598775598299\" "
                        "theta=\"0.785398163397448\" />");
        Angle a(xml.child("angle"));

        CHECK_CLOSE(0.523598775598299, a.alpha, 0.000000000001);
        CHECK_CLOSE(0.785398163397448, a.theta, 0.000000000001);
        CHECK_CLOSE(0.612372435695794, a.ox, 0.00000000001);
        CHECK_CLOSE(0.353553390593274, a.oy, 0.00000000001);
        CHECK_CLOSE(0.707106781186548, a.oz, 0.00000000001);
        CHECK_EQUAL(0.1, a.weight);
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
