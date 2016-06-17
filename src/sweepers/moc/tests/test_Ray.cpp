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

#include <cassert>
#include <string>
#include "pugixml.hpp"
#include "util/global_config.hpp"
#include "core/angular_quadrature.hpp"
#include "core/constants.hpp"
#include "core/core_mesh.hpp"
#include "moc/ray_data.hpp"

using namespace mocc;
using namespace mocc::moc;

using std::cout;
using std::endl;

TEST(simple_ray)
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("6x5.xml");

    CHECK(result);

    mocc::CoreMesh mesh(geom_xml);
    {

        // Test a few rays that starts on a corner, ends on a corner and crosses
        // a bunch of corners
        {
            Ray ray(Point2(0.0, 1.0), Point2(4.0, 5.0), {0, 0}, 0, mesh);

            CHECK_EQUAL(ray.cm_surf_fw(), 37);
            CHECK_EQUAL(ray.cm_cell_fw(), 6);
            CHECK_EQUAL(ray.cm_surf_bw(), 88);
            CHECK_EQUAL(ray.cm_cell_bw(), 27);

            CHECK_EQUAL(ray.nseg(), 12);
            CHECK_EQUAL(ray.ncseg(), 8);

            // all of the segment lengths should be the same. I'm not testing
            // this too much in the general sense, since the tests for the pin
            // meshes should find most of these types of issues.
            real_t t = 1.0 / 3.0 * sqrt(2);
            for (auto v : ray.seg_len()) {
                CHECK_CLOSE(v, t, 0.00001);
            }

            std::vector<Surface> fw_surf = {
                Surface::EAST, Surface::NORTH, Surface::EAST, Surface::NORTH,
                Surface::EAST, Surface::NORTH, Surface::EAST, Surface::NORTH};

            std::vector<Surface> bw_surf = {
                Surface::WEST, Surface::SOUTH, Surface::WEST,  Surface::SOUTH,
                Surface::WEST, Surface::SOUTH, Surface::SOUTH, Surface::WEST};

            VecI nseg = {3, 0, 3, 0, 3, 0, 3, 0};
            for (int i = 0; i < ray.ncseg(); i++) {
                auto rcd = ray.cm_data()[i];
                CHECK_EQUAL(rcd.fw, fw_surf[i]);
                CHECK_EQUAL(rcd.bw, bw_surf[i]);
                CHECK_EQUAL(rcd.nseg_fw, nseg[i]);
                CHECK_EQUAL(rcd.nseg_bw, nseg[i]);
            }
        }

        {
            Ray ray(Point2(4.0, 0.0), Point2(6.0, 2.0), {0, 0}, 0, mesh);

            CHECK_EQUAL(ray.cm_surf_fw(), 89);
            CHECK_EQUAL(ray.cm_cell_fw(), 4);
            CHECK_EQUAL(ray.cm_surf_bw(), 43);
            CHECK_EQUAL(ray.cm_cell_bw(), 11);

            CHECK_EQUAL(ray.nseg(), 6);
            CHECK_EQUAL(ray.ncseg(), 4);
        }

        {
            Ray ray(Point2(2.0, 0.0), Point2(0.0, 2.0), {0, 0}, 0, mesh);
        }

        {
            Ray ray(Point2(6.0, 3.0), Point2(4.0, 5.0), {0, 0}, 0, mesh);
        }

        {
            Ray ray(Point2(0.0, 0.5), Point2(6.0, 3.25), {0, 0}, 0, mesh);
        }
    }
}

TEST(nasty_ray)
{

    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("square.xml");

    mocc::CoreMesh mesh(geom_xml);

    pugi::xml_document angquad_xml;
    result = angquad_xml.load_string("<ang_quad type=\"ls\" order=\"4\" />");
    // Make a nasty ray to exercise the coarse indexing
    {
        Ray ray(Point2(1.26, 0.0), Point2(3.78, 2.52), {0, 0}, 0, mesh);
    }
    {
        Ray ray(Point2(1.26, 0.0), Point2(0.0, 1.26), {0, 0}, 0, mesh);
    }
    {
        Ray ray(Point2(0.0, 1.26), Point2(2.52, 3.78), {0, 0}, 0, mesh);
    }
    {
        Ray ray(Point2(3.78, 2.52), Point2(2.52, 3.78), {0, 0}, 0, mesh);
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
