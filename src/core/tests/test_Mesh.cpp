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
#include <vector>

#include "core/constants.hpp"
#include "core/fp_utils.hpp"
#include "core/geom.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"

using std::cout;
using std::endl;

using namespace mocc;

TEST( mesh )
{
    // Make a simple mesh, 1.0 cm pitch to keep things simple, 6x5
    mocc::VecF x;
    mocc::VecF y;
    mocc::VecF z;
    for( float xi=0.0; xi<6.1; xi++ ) {
        x.push_back(xi);
    }
    for( float yi=0.0; yi<5.1; yi++ ) {
        y.push_back(yi);
    }

    z.push_back(0.0);
    z.push_back(1.0);
    z.push_back(3.0);

    std::array<Boundary, 6> bc =
    {
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT
    };

    mocc::Mesh mesh( 30, 30, x, y, z, bc );


    // Test the boundary point stuff (coarse_boundary_cell())
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(0.0, 2.0), 1), 12 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(0.0, 3.0), 4), 12 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(2.0, 5.0), 4), 26 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(4.0, 5.0), 3), 27 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(6.0, 4.0), 3), 23 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(6.0, 2.0), 2), 17 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(2.0, 0.0), 2), 1 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(3.5, 0.0), 1), 3 );
    CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(5.0, 0.0), 1), 5 );


    // Test surface normals
    CHECK_EQUAL( mesh.surface_normal(30), Normal::X_NORM );
    CHECK_EQUAL( mesh.surface_normal(47), Normal::X_NORM );
    CHECK_EQUAL( mesh.surface_normal(57), Normal::X_NORM );
    CHECK_EQUAL( mesh.surface_normal(64), Normal::X_NORM );
    CHECK_EQUAL( mesh.surface_normal(58), Normal::X_NORM );

    CHECK_EQUAL( mesh.surface_normal(69), Normal::Y_NORM );
    CHECK_EQUAL( mesh.surface_normal(100), Normal::Y_NORM );
    CHECK_EQUAL( mesh.surface_normal(65), Normal::Y_NORM );
    CHECK_EQUAL( mesh.surface_normal(95), Normal::Y_NORM );
    CHECK_EQUAL( mesh.surface_normal(74), Normal::Y_NORM );
    CHECK_EQUAL( mesh.surface_normal(70), Normal::Y_NORM );

    CHECK_EQUAL( mesh.surface_normal(0), Normal::Z_NORM );
    CHECK_EQUAL( mesh.surface_normal(29), Normal::Z_NORM );
    CHECK_EQUAL( mesh.surface_normal(14), Normal::Z_NORM );
    CHECK_EQUAL( mesh.surface_normal(101), Normal::Z_NORM );
    CHECK_EQUAL( mesh.surface_normal(129), Normal::Z_NORM );

    // Test cells straddling surfaces
    // X normals
    CHECK_EQUAL( mesh.coarse_neigh_cells(53).first,  19);
    CHECK_EQUAL( mesh.coarse_neigh_cells(53).second, 20);
    CHECK_EQUAL( mesh.coarse_neigh_cells(37).first, -1);
    CHECK_EQUAL( mesh.coarse_neigh_cells(37).second, 6);
    CHECK_EQUAL( mesh.coarse_neigh_cells(64).first, 29);
    CHECK_EQUAL( mesh.coarse_neigh_cells(64).second, -1);

    // Y normals
    CHECK_EQUAL( mesh.coarse_neigh_cells(65).first, -1 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(65).second, 0);
    CHECK_EQUAL( mesh.coarse_neigh_cells(100).first, 29 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(100).second, -1);
    CHECK_EQUAL( mesh.coarse_neigh_cells(80).first, 14 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(80).second, 20);
    CHECK_EQUAL( mesh.coarse_neigh_cells(83).first, -1 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(83).second, 3);

    // Z normals
    CHECK_EQUAL( mesh.coarse_neigh_cells(0).first, -1 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(0).second, 0 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(29).first, -1 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(29).second, 29 );
    CHECK_EQUAL( mesh.coarse_neigh_cells(115).first, 14);
    CHECK_EQUAL( mesh.coarse_neigh_cells(115).second, 44 );
}


// Test a more irregular mesh. make sure the volume and area stuff comes out
// okay.
TEST( test_irregular )
{
    // Make a simple mesh, 1.0 cm pitch to keep things simple, 6x5
    mocc::VecF x = {0.0, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0};
    mocc::VecF y = {0.0, 1.0, 2.0, 3.5, 4.0, 4.5, 7.0};
    mocc::VecF z;

    z.push_back(0.0);
    z.push_back(1.0);
    z.push_back(3.0);

    std::array<Boundary, 6> bc =
    {
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT
    };

    mocc::Mesh mesh( 30, 30, x, y, z, bc );

    CHECK_EQUAL( mesh.coarse_area(78), 1.0 );
    CHECK_EQUAL( mesh.coarse_area(83), 1.0 );
    CHECK_EQUAL( mesh.coarse_area(87), 1.0 );
    CHECK_EQUAL( mesh.coarse_area(91), 1.0 );
    CHECK_EQUAL( mesh.coarse_area(93), 0.5 );
    CHECK_EQUAL( mesh.coarse_area(95), 0.5 );
    CHECK_EQUAL( mesh.coarse_area(71), 2.5 );
    CHECK_EQUAL( mesh.coarse_area(77), 2.5 );
    CHECK_EQUAL( mesh.coarse_area(64), 0.5 );
    CHECK_EQUAL( mesh.coarse_area(60), 0.5 );
    CHECK_EQUAL( mesh.coarse_area(14), 0.75 );
    CHECK_EQUAL( mesh.coarse_area(31), 2.5 );
    CHECK_EQUAL( mesh.coarse_area(32), 1.25 );
}

int main() {
    return UnitTest::RunAllTests();
}
