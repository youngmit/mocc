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

#include "pugixml.hpp"
#include "core/assembly.hpp"
#include "core/lattice.hpp"

#include "inputs.hpp"

using namespace mocc;

TEST(assembly)
{
    std::string test_xml = "<mesh type=\"rect\" id=\"1\" pitch=\"1.2\">"
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
                           "<lattice id=\"2\" nx=\"3\" ny=\"5\">"
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
                           "       2"
                           "       1"
                           "       2"
                           "       2"
                           "   </lattices>"
                           "</assembly>"
                           "<assembly id=\"2\" np=\"5\">"
                           "   <hz>"
                           "       3.14 3.14 3.14 3.14 3.14"
                           "   </hz>"
                           "   <lattices>"
                           "       1"
                           "       2"
                           "       1"
                           "       2"
                           "       2"
                           "   </lattices>"
                           "</assembly>"
                           ""
                           "<assembly id=\"3\" np=\"5\">"
                           "   <hz>"
                           "       1.0 1.0 1.0 2.0 1.0"
                           "   </hz>"
                           "   <lattices>"
                           "       1"
                           "       1"
                           "       1"
                           "       1"
                           "       1"
                           "   </lattices>"
                           "</assembly>"
                           "<assembly id=\"4\" np=\"20\" hz=\"1.23\">"
                           "   <lattices>"
                           "      { 1 1 1 1 1 } "
                           "      { 2 2 2 2 2 "
                           "       2 2 2 2 2 }"
                           "       1 1 2 2 2"
                           "   </lattices>"
                           "</assembly>"
                           "";

    pugi::xml_document xml;
    auto result = xml.load_string(test_xml.c_str());

    REQUIRE CHECK(result);

    auto meshes = ParsePinMeshes(xml);
    MaterialLib mat_lib(xml.child("material_lib"));
    auto pins       = ParsePins(xml, meshes, mat_lib);
    auto lattices   = ParseLattices(xml, pins);
    auto assemblies = ParseAssemblies(xml, lattices);

    CHECK(assemblies[1]->compatible(*assemblies[2]));
    CHECK(!assemblies[1]->compatible(*assemblies[3]));

    std::vector<int> ref_subplane = {1, 1, 1, 1, 1, 10, 5};
    CHECK_EQUAL(7, assemblies[4]->subplane().size());
    CHECK_ARRAY_EQUAL(ref_subplane, assemblies[4]->subplane(),7);
}

TEST(more) {
    std::string xml =
        pinmesh_xml + material_xml + pin_xml + lattice_xml + assembly_xml;

    pugi::xml_document xml_doc;
    pugi::xml_parse_result result = xml_doc.load_string(xml.c_str());

    REQUIRE CHECK(result);

    auto meshes = ParsePinMeshes(xml_doc);
    MaterialLib mat_lib(xml_doc.child("material_lib"));
    auto pins = ParsePins(xml_doc, meshes, mat_lib);
    auto lattices = ParseLattices(xml_doc, pins);
    auto assemblies = ParseAssemblies(xml_doc, lattices);

    CHECK(assemblies[1]->compatible(*assemblies[2]));
    CHECK(!assemblies[1]->compatible(*assemblies[100]));

}

int main(int, const char *[])
{
    return UnitTest::RunAllTests();
}
