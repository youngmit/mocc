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

#include "xs_mesh_homogenized.hpp"

#include "inputs.hpp"

using namespace mocc;
using std::cout;
using std::endl;

std::string core_xml = "<core nx=\"4\" ny=\"3\""
                       "    north  = \"reflect\""
                       "    south  = \"reflect\""
                       "    east   = \"reflect\""
                       "    west   = \"reflect\""
                       "    top    = \"vacuum\""
                       "    bottom = \"vacuum\" > "
                       "    1 1 1 1"
                       "    1 1 1 1"
                       "    1 1 1 1"
                       "</core>";

TEST(xsmeshhom)
{
    {
        pugi::xml_document geom_xml;
        pugi::xml_parse_result result = geom_xml.load_file("2x3_1.xml");
        CHECK(result);

        CoreMesh mesh(geom_xml);

        XSMeshHomogenized xs_mesh(mesh);

        H5Node h5f("xsmesh_1.h5", H5Access::WRITE);
        xs_mesh.output(h5f);
    }
    {
        pugi::xml_document geom_xml;
        pugi::xml_parse_result result = geom_xml.load_file("2x3_2.xml");
        CHECK(result);

        CoreMesh mesh(geom_xml);

        XSMeshHomogenized xs_mesh(mesh);

        H5Node h5f("xsmesh_2.h5", H5Access::WRITE);
        xs_mesh.output(h5f);
    }
}

// Tests some of the error checking involved in constructing an XSMeshHom from
// data files.
TEST(fromdata_fail)
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("2x3_stack.xml");
    CHECK(result);
    if (!result) {
        std::cout << result.description() << std::endl;
    }
    // the above xml file doesnt come with a core tag, so add one
    // assembly 1 is a regular stack of lattices
    result = geom_xml.append_buffer(core_xml.c_str(), core_xml.size());

    if (!result) {
        std::cout << "error parsing <core> addendum: " << std::endl;
        std::cout << result.description() << std::endl;
    }

    CoreMesh mesh(geom_xml);

    // Test error checks on the XML input
    {
        cout << "invalid" << endl;
        pugi::xml_document xsmesh_xml;
        std::string xml = "<data file=\"xsmesh_1.h5\" top_plane=\"-1\"/>";
        xsmesh_xml.load_string(xml.c_str());

        CHECK_THROW(XSMeshHomogenized xs_mesh(mesh, xsmesh_xml), Exception);
    }
    {
        cout << "out of order" << endl;
        pugi::xml_document xsmesh_xml;
        std::string xml = "<data file=\"xsmesh_1.h5\" top_plane=\"5\"/>"
                          "<data file=\"xsmesh_2.h5\" top_plane=\"1\"/>";
        xsmesh_xml.load_string(xml.c_str());

        CHECK_THROW(XSMeshHomogenized xs_mesh(mesh, xsmesh_xml), Exception);
    }
}

// Test an actual xsmesh hom object that should successfully construct
TEST(fromdata)
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("2x3_stack.xml");
    CHECK(result);
    if (!result) {
        std::cout << result.description() << std::endl;
    }
    result = geom_xml.append_buffer(core_xml.c_str(), core_xml.size());
    if (!result) {
        std::cout << "error parsing <core> addendum: " << std::endl;
        std::cout << result.description() << std::endl;
    }

    CoreMesh mesh(geom_xml);

    pugi::xml_document xsmesh_xml;
    std::string xml =
        "<data file=\"xsmesh_2.h5\" bottom_plane=\"0\" top_plane=\"7\"/>"
        "<data file=\"xsmesh_1.h5\" bottom_plane=\"8\" top_plane=\"11\"/>"
        "";
    xsmesh_xml.load_string(xml.c_str());

    XSMeshHomogenized xs_mesh(mesh, xsmesh_xml);

    CHECK_EQUAL(7, xs_mesh.eubounds().size());

    CHECK_EQUAL(864, xs_mesh.size());

    for(int i=576; i<864; i++) {
        CHECK_CLOSE(2.005998E-02, xs_mesh[i].xsmacnf(0), 0.000001);
    }

    for(int ilat=0; ilat<96; ilat++) {
        CHECK_CLOSE(0.0125521, xs_mesh[ilat*6+0].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0125521, xs_mesh[ilat*6+1].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0115752, xs_mesh[ilat*6+2].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0115752, xs_mesh[ilat*6+3].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0115752, xs_mesh[ilat*6+4].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0125521, xs_mesh[ilat*6+5].xsmacnf(0), 0.000001);
    }

    // test this whole plane. It's got some asymmetry, so if everything is good
    // here, we can be pretty certain of the X- and Y-dimensions in the transfer
    CHECK_CLOSE(0.0125521, xs_mesh[36].xsmacnf(0), 0.000001);
    CHECK_CLOSE(0.0125521, xs_mesh[37].xsmacnf(0), 0.000001);
    CHECK_CLOSE(0.0115752, xs_mesh[38].xsmacnf(0), 0.000001);
    CHECK_CLOSE(0.0115752, xs_mesh[39].xsmacnf(0), 0.000001);
    CHECK_CLOSE(0.0115752, xs_mesh[40].xsmacnf(0), 0.000001);
    CHECK_CLOSE(0.0125521, xs_mesh[41].xsmacnf(0), 0.000001);

    // Now for the big guns: Re-make the cross section mesh directly from the
    // core mesh via homogenization and check all of the fields. This doesnt
    // test the actual homogenization procedures, but its an excellent test of
    // the I/O procedures.
    XSMeshHomogenized xs_reference(mesh);

    CHECK(xs_mesh == xs_reference);
}

// Test creation of an XS mesh using the macroplane grouping. Use assembly 2 in
// 2x3_stack.
TEST(macroplanes)
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("2x3_stack.xml");
    REQUIRE CHECK(result);
    if (!result) {
        std::cout << result.description() << std::endl;
    }
    // the above xml file doesnt come with a core tag, so add one
    // assembly 1 is a regular stack of lattices
    std::string core_xml = "<core nx=\"4\" ny=\"3\""
                           "    north  = \"reflect\""
                           "    south  = \"reflect\""
                           "    east   = \"reflect\""
                           "    west   = \"reflect\""
                           "    top    = \"vacuum\""
                           "    bottom = \"vacuum\" > "
                           "    2 2 2 2 "
                           "    2 2 2 2 "
                           "    2 2 2 2 "
                           "</core>";
    result = geom_xml.append_buffer(core_xml.c_str(), core_xml.size());
    if (!result) {
        std::cout << "error parsing <core> addendum: " << std::endl;
        std::cout << result.description() << std::endl;
    }

    CoreMesh mesh(geom_xml);

    XSMeshHomogenized xs_mesh(mesh);

    CHECK_EQUAL(288, xs_mesh.size());

    for(int i=216; i<288; i++) {
        CHECK_CLOSE(2.005998E-02, xs_mesh[i].xsmacnf(0), 0.000001);
    }
    for(int ilat=0; ilat<36; ilat++) {
        CHECK_CLOSE(0.0125521, xs_mesh[ilat*6+0].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0125521, xs_mesh[ilat*6+1].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0115752, xs_mesh[ilat*6+2].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0115752, xs_mesh[ilat*6+3].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0115752, xs_mesh[ilat*6+4].xsmacnf(0), 0.000001);
        CHECK_CLOSE(0.0125521, xs_mesh[ilat*6+5].xsmacnf(0), 0.000001);
    }
}

TEST(complicated) {
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_string(complex_xml.c_str());
    REQUIRE CHECK(result);

    CoreMesh mesh(geom_xml);
    XSMeshHomogenized xs_mesh(mesh);

    CHECK_EQUAL(960, xs_mesh.size());
    
    VecF xs(mesh.n_reg(MeshTreatment::PIN), 0.0);
    int ixs = 0;
    for(const auto &xsreg: xs_mesh) {
        real_t xs_i = xsreg.xsmactr(0);
        std::cout << ixs << " : ";
        for(const int reg: xsreg.reg()){
            std::cout << reg << " ";
            xs[reg] = xs_i;
        }
        ixs++;
        std::cout << std::endl;
    }

    int icell = mesh.coarse_cell(Position(1,6,19));
    std::cout << icell << std::endl;
    CHECK_CLOSE(1.40063733419359E-01, xs[icell], REAL_FUZZ);
    icell = mesh.coarse_cell(Position(1,6,20));
    CHECK_CLOSE(1.9242061690222E-01, xs[icell], REAL_FUZZ);
}

int main()
{
    return UnitTest::RunAllTests();
}
