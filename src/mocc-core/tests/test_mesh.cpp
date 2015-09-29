#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <cassert>
#include <iostream>
#include <vector>

#include "constants.hpp"
#include "fp_utils.hpp"
#include "geom.hpp"
#include "global_config.hpp"
#include "mesh.hpp"
#include "ray_data.hpp"

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
    for( float xi=0.0; xi<6.1; xi++ ) {
        x.push_back(xi);
    }
    for( float yi=0.0; yi<5.1; yi++ ) {
        y.push_back(yi);
    }

    mocc::Mesh mesh( 30, 30, 6, 5, 1, x, y);


    // Test the boundary point stuff (coarse_boundary_cell())
    cout << "hi there" << endl;
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(0.0, 2.0), 1) == 12 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(0.0, 3.0), 4) == 12 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(2.0, 5.0), 4) == 26 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(4.0, 5.0), 3) == 27 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(6.0, 4.0), 3) == 23 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(6.0, 2.0), 2) == 17 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(2.0, 0.0), 2) == 1 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(3.5, 0.0), 1) == 3 );
    BOOST_CHECK( mesh.coarse_boundary_cell(Point2(5.0, 0.0), 1) == 5 );


    /*
     * disabled since the conventions for how coarse mesh intersections are
     * performed have changed.
    // check surface indexing
    int s[2];
    Point2 p(0.5, 0.0);
    Point2 pcell(0.5, 0.5);
    int cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK( cell == 0 );
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 65);

    // check a boundary corner on the south face
    p.x = 3.0; p.y = 0.0;
    pcell.x = 3.5; pcell.y = 0.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 3);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 83 );
    // boundary corner on north face
    p.x = 4.0; p.y = 5.0;
    pcell.x = 4.5; pcell.y = 4.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 28);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 94 );
    // boundary corner on east face
    p.x = 6.0; p.y = 4.0;
    pcell.x = 5.5; pcell.y = 3.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 23);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 57 );
    // boundary corner on west face
    p.x = 0.0; p.y = 2.0;
    pcell.x = 0.5; pcell.y = 2.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 12);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 44 );
    // interior x normal face
    p.x = 3.0; p.y = 1.5;
    pcell.x = 2.5; pcell.y = 1.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 8);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 40 );
    // interior y normal face
    p.x = 2.5; p.y = 4.0;
    pcell.x = 2.5; pcell.y = 3.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 20);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 1 );
    BOOST_CHECK(s[0] == 81 );
    // interior corner
    p.x = 4.0; p.y = 2.0;
    pcell.x = 4.5; pcell.y = 2.5;
    cell = mesh.coarse_cell_point( pcell );
    BOOST_CHECK(cell == 16);
    BOOST_CHECK(mesh.coarse_surf_point( p, cell, s ) == 2 );
    BOOST_CHECK(s[0] == 48 );
    BOOST_CHECK(s[1] == 85 );

    */

}
