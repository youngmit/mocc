#undef __STRICT_ANSI__
#undef _REENT_ONLY
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <cassert>
#include <iostream>
#include <vector>

#include "mocc-core/constants.hpp"
#include "mocc-core/fp_utils.hpp"
#include "mocc-core/geom.hpp"
#include "mocc-core/global_config.hpp"
#include "mocc-core/mesh.hpp"

using std::cout;
using std::endl;

using namespace mocc;

BOOST_AUTO_TEST_CASE( testall )
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

    Boundary bc[6] =
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
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(0.0, 2.0), 1), 12 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(0.0, 3.0), 4), 12 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(2.0, 5.0), 4), 26 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(4.0, 5.0), 3), 27 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(6.0, 4.0), 3), 23 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(6.0, 2.0), 2), 17 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(2.0, 0.0), 2), 1 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(3.5, 0.0), 1), 3 );
    BOOST_CHECK_EQUAL( mesh.coarse_boundary_cell(Point2(5.0, 0.0), 1), 5 );


    // Test surface normals
    BOOST_CHECK_EQUAL( mesh.surface_normal(30), Normal::X_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(47), Normal::X_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(57), Normal::X_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(64), Normal::X_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(58), Normal::X_NORM );

    BOOST_CHECK_EQUAL( mesh.surface_normal(69), Normal::Y_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(100), Normal::Y_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(65), Normal::Y_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(95), Normal::Y_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(74), Normal::Y_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(70), Normal::Y_NORM );

    BOOST_CHECK_EQUAL( mesh.surface_normal(0), Normal::Z_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(29), Normal::Z_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(14), Normal::Z_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(101), Normal::Z_NORM );
    BOOST_CHECK_EQUAL( mesh.surface_normal(129), Normal::Z_NORM );

    // Test cells straddling surfaces
    // X normals
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(53).first,  19);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(53).second, 20);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(37).first, -1);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(37).second, 6);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(64).first, 29);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(64).second, -1);

    // Y normals
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(65).first, -1 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(65).second, 0);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(100).first, 29 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(100).second, -1);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(80).first, 14 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(80).second, 20);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(83).first, -1 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(83).second, 3);

    // Z normals
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(0).first, -1 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(0).second, 0 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(29).first, -1 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(29).second, 29 );
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(115).first, 14);
    BOOST_CHECK_EQUAL( mesh.coarse_neigh_cells(115).second, 44 );


}


// Test a more irregular mesh. make sure the volume and area stuff comes out
// okay.
BOOST_AUTO_TEST_CASE( test_irregular )
{
    // Make a simple mesh, 1.0 cm pitch to keep things simple, 6x5
    mocc::VecF x = {0.0, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0};
    mocc::VecF y = {0.0, 1.0, 2.0, 3.5, 4.0, 4.5, 7.0};
    mocc::VecF z;

    z.push_back(0.0);
    z.push_back(1.0);
    z.push_back(3.0);

    Boundary bc[6] =
    {
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT
    };

    mocc::Mesh mesh( 30, 30, x, y, z, bc );

    BOOST_CHECK_EQUAL( mesh.coarse_area(78), 1.0 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(83), 1.0 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(87), 1.0 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(91), 1.0 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(93), 0.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(95), 0.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(71), 2.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(77), 2.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(64), 0.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(60), 0.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(14), 0.75 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(31), 2.5 );
    BOOST_CHECK_EQUAL( mesh.coarse_area(32), 1.25 );


}
