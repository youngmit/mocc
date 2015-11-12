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
#include "mocc-core/ray_data.hpp"

using std::cout;
using std::endl;

using namespace mocc;

BOOST_AUTO_TEST_CASE( testall )
{
    // For now, mostly just testing the coarse ray data, since its the easiest
    // to screw up.

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

}
