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

#include "core/assembly.hpp"
#include "core/lattice.hpp"


TEST( assembly )
{
    std::string test_xml =
        "<mesh type=\"rec\" pitch=\"1.2\">"
        "<sub_x>1</sub_x>"
        "<sub_y>1</sub_y>"
        "</mesh>"
        ""
        "<material_lib path=\"c5g7.xsl\">"
        "   <material id=\"1\" name=\"UO2-3.3\" />"
        "</material_lib>"
        ""
        "<pin id=\"1\" mesh=\"1\">"
        "   1"
        "</pin>"
        ""
        "<lattice id=\"1\" nx=\"3\" ny=\"5\">"
        "   1 1 1"
        "   1 1 1"
        "   1 1 1"
        "   1 1 1"
        "   1 1 1"
        "</lattice>"
        ""
        "<assembly id=\"1\" np=\"5\" hz=\"3.14\">"
        "   <lattices>"
        "       1"
        "       1"
        "       1"
        "       1"
        "       1"
        "   </lattices>"
        "</assembly>"
        ""
        "<assembly id=\"1\" np=\"5\" hz=\"3.14\">"
        "   <lattices>"
        "       1"
        "       1"
        "       1"
        "       1"
        "       1"
        "   </lattices>"
        "</assembly>"
        "";

    pugi::xml_document xml;
    auto result = xml.load_string( test_xml.c_str() );

    CHECK(result);


}

int main(int, const char*[]) {
    return UnitTest::RunAllTests();
}
