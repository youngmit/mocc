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
#include "util/h5file.hpp"
#include "core/core_mesh.hpp"
#include "sweepers/mc/particle_pusher.hpp"
#include "sweepers/mc/particle.hpp"
#include "util/rng_lcg.hpp"

using namespace mocc;
using namespace mocc::mc;

TEST(test_point)
{
    pugi::xml_document geom_xml;
    geom_xml.load_file("square.xml");

    CoreMesh mesh(geom_xml);

    XSMesh xs_mesh(mesh);

    ParticlePusher pusher(mesh, xs_mesh);

    int N = 10000;
    int pc = N/100;
    RNG_LCG rng(11112854149);
    for (int i = 0; i < N; i++) {
        if(i % pc == 0) {
            std::cout << i/pc << "%" << std::endl;
        }
        Point3 loc(0.00000001, 0.000000001, 0.25);
        Direction dir(rng.random(TWOPI), rng.random(PI));
        Particle p(loc, dir, 0, i);

        pusher.simulate(p, true);
    }

    H5Node h5("point_source.h5", H5Access::WRITE);
    pusher.output(h5);

    return;
}

TEST(test_beam)
{
    pugi::xml_document geom_xml;
    geom_xml.load_file("tunnel.xml");

    CoreMesh mesh(geom_xml);

    XSMesh xs_mesh(mesh);

    ParticlePusher pusher(mesh, xs_mesh);

    int N = 10000;
    int pc = N/100;
    RNG_LCG rng(11112854149);
    for (int i = 0; i < N; i++) {
        if(i % pc == 0) {
            std::cout << i/pc << "%" << std::endl;
        }
        Point3 loc(0.00000000000001, rng.random(4.5), rng.random(0.5));
        Direction dir(0.0, HPI);
        Particle p(loc, dir, 0, i);

        pusher.simulate(p, true);
    }

    H5Node h5("exp.h5", H5Access::WRITE);
    pusher.output(h5);

    return;
}

int main()
{
    return UnitTest::RunAllTests();
}
