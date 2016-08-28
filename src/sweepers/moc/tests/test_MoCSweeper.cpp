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
#include <string>
#include "pugixml.hpp"
#include "util/blitz_typedefs.hpp"
#include "util/global_config.hpp"
#include "core/core_mesh.hpp"
#include "core/eigen_interface.hpp"
#include "core/material_lib.hpp"
#include "sweepers/moc/moc_sweeper.hpp"

#include "core/tests/inputs.hpp"

using namespace mocc;
using moc::MoCSweeper;

TEST(moc)
{
    pugi::xml_document xml_doc;
    {
        auto result = xml_doc.load_string(complex_xml.c_str());
        CHECK(result);
    }

    CoreMesh mesh(xml_doc);

    pugi::xml_document moc_doc;
    {
        std::string moc_xml = "<sweeper type=\"moc\" n_inner=\"10\">"
                              "    <ang_quad type=\"ls\" order=\"6\"/>"
                              "    <rays spacing=\"0.05\"/>"
                              "</sweeper>";
        auto result = moc_doc.load_string(moc_xml.c_str());
        CHECK(result);
    }

    MoCSweeper sweeper(moc_doc.child("sweeper"), mesh);

    sweeper.initialize();

    // Make sure that we can get/set/get pin flux without change
    {
        ArrayB1 pin_flux_1(mesh.n_reg(MeshTreatment::PLANE));
        ArrayB1 pin_flux_2(mesh.n_reg(MeshTreatment::PLANE));

        sweeper.get_pin_flux_1g(0, pin_flux_1);
        sweeper.set_pin_flux_1g(0, pin_flux_1);
        sweeper.get_pin_flux_1g(0, pin_flux_2);

        bool good = true;
        for (int i = 0; i < pin_flux_1.size(); i++) {
            CHECK_CLOSE(pin_flux_1(i), pin_flux_2(i), REAL_FUZZ);
            if (!fp_equiv_ulp(pin_flux_1(i), pin_flux_2(i))) {
                good = false;
                std::cout << i << " " << pin_flux_1(i) << " " << pin_flux_2(i)
                          << "\n";
            }
        }
    }
    {
        ArrayB1 pin_flux_1(mesh.n_reg(MeshTreatment::PIN));
        ArrayB1 pin_flux_2(mesh.n_reg(MeshTreatment::PIN));

        sweeper.get_pin_flux_1g(0, pin_flux_1, MeshTreatment::PIN);
        sweeper.set_pin_flux_1g(0, pin_flux_1, MeshTreatment::PIN);
        sweeper.get_pin_flux_1g(0, pin_flux_2, MeshTreatment::PIN);

        bool good = true;
        for (int i = 0; i < pin_flux_1.size(); i++) {
            CHECK_CLOSE(pin_flux_1(i), pin_flux_2(i), REAL_FUZZ);
            if (!fp_equiv_ulp(pin_flux_1(i), pin_flux_2(i))) {
                good = false;
                std::cout << i << " " << pin_flux_1(i) << " " << pin_flux_2(i)
                          << "\n";
            }
        }
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
