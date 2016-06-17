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

#include <memory>

#include "pugixml.hpp"

#include "core/tests/pugi_utils.hpp"

#include "core/cmfd.hpp"
#include "core/xs_mesh_homogenized.hpp"

using namespace mocc;

/// \todo this test is practically nonexistent

TEST( testCMFD ) {
    auto mesh_xml = inline_xml_file("3x5.xml");
    CoreMesh mesh(*mesh_xml);

    auto cmfd_xml = inline_xml("<cmfd k_tol=\"1e-10\" "
            "psi_tol=\"1e-8\" "
            "max_iter=\"100\" "
            "enabled=\"t\" "
            "negative_fixup=\"f\" />");

    std::shared_ptr<XSMeshHomogenized> xsmesh(
            std::make_shared<XSMeshHomogenized>( mesh ) );

    CMFD cmfd(*cmfd_xml, &mesh, xsmesh);

    real_t k = 1.0;
    cmfd.solve( k );
    std::cout << k << std::endl;
}


int main() {
    return UnitTest::RunAllTests();
}
