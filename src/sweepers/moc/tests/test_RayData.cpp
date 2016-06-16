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
#include <iostream>
#include <string>
#include "pugixml.hpp"
#include "util/global_config.hpp"
#include "angular_quadrature.hpp"
#include "constants.hpp"
#include "core_mesh.hpp"
#include "ray_data.hpp"

using namespace mocc;

TEST(raydata)
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("square.xml");

    CoreMesh mesh(geom_xml);
    std::cout << "mesh" << std::endl;

    pugi::xml_document angquad_xml;
    result = angquad_xml.load_string("<ang_quad type=\"ls\" order=\"4\" />");

    CHECK(result);

    AngularQuadrature ang_quad(angquad_xml.child("ang_quad"));

    pugi::xml_document ray_xml;
    ray_xml.load_string("<rays spacing=\"0.01\" />");

    moc::RayData ray_data(ray_xml.child("rays"), ang_quad, mesh);

    for (auto &plane_rays : ray_data) {
        int iang    = 0;
        double wsum = 0.0;
        VecF vol(mesh.n_reg(), 0.0);
        for (auto &angle_rays : plane_rays) {

            wsum += ang_quad[iang].weight * 2.0 * PI;
            real_t space  = ray_data.spacing(iang);
            real_t wt_ang = space * ang_quad[iang].weight * 2.0 * PI;
            for (auto &ray : angle_rays) {
                for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                    int ireg = ray.seg_index(iseg);
                    vol[ireg] += ray.seg_len(iseg) * wt_ang;
                }
            }

            iang++;
        }
        for (auto &v : vol) {
            CHECK_CLOSE(0.1764 * 4.0 * PI, v, 0.00000000000001);
        }
        std::cout << "angle weight sum: " << wsum << std::endl;
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
