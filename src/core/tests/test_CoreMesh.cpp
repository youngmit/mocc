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
#include "core/core_mesh.hpp"

#include "inputs.hpp"

using namespace mocc;

TEST(assembly)
{
    pugi::xml_document xml_doc;
    pugi::xml_parse_result result = xml_doc.load_string(complex_xml.c_str());

    REQUIRE CHECK(result);

    CoreMesh mesh(xml_doc);

    {
        const MacroPlane &macroplane = mesh.macroplanes()[3];
        const Pin &pin = **(macroplane.begin()+44);
        CHECK(pin == *(mesh.pins().at(7)));
    }
    {
        const MacroPlane &macroplane = mesh.macroplanes()[4];
        const Pin &pin = **(macroplane.begin()+44);
        CHECK(pin == *(mesh.pins().at(8)));
    }
    {
        const MacroPlane &macroplane = mesh.macroplanes()[6];
        const Pin &pin = **(macroplane.begin()+44);
        CHECK(pin == *(mesh.pins().at(8)));
    }
    {
        const MacroPlane &macroplane = mesh.macroplanes()[7];
        const Pin &pin = **(macroplane.begin()+45);
        CHECK(pin == *(mesh.pins().at(1)));
    }
    {
        const MacroPlane &macroplane = mesh.macroplanes()[8];
        const Pin &pin = **(macroplane.begin()+45);
        CHECK(pin == *(mesh.pins().at(6)));
    }

    for(const auto mp: mesh.macroplanes()) {
        std::cout << mp << std::endl;
    }

    std::cout << mesh.macroplanes().front().back()->id() << std::endl;

}

int main()
{
    UnitTest::RunAllTests();
}
