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

#include "xs_mesh.hpp"

using namespace std;
using namespace mocc;

TEST(test)
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file("2x3_1.xml");
    CHECK(result);

    CoreMesh mesh(geom_xml);

    XSMesh xs_mesh(mesh);

    cout << "Sig-t: " << xs_mesh[0].xsmactr(0) << endl;

    auto cdf = xs_mesh[0].reaction_cdf(0);
    cout << "Reaction CDF:" << endl;
    for (const auto &v : cdf) {
        cout << v << endl;
    }

    cdf = xs_mesh[0].chi_cdf();
    cout << "Chi CDF:" << endl;
    for (const auto &v : cdf) {
        cout << v << endl;
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
