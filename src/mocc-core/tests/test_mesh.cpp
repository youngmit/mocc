#include <iostream>
#include <cassert>
#include <vector>

#include "global_config.hpp"
#include "geom.hpp"
#include "constants.hpp"
#include "mesh.hpp"
#include "ray_data.hpp"
#include "fp_utils.hpp"

using std::cout;
using std::endl;

using mocc::Point2;

int main() {
    // For now, mostly just testing the coarse ray data, since its the easiest
    // to screw up.
    
    // Make a simple mesh, 1.0 cm pitch to keep things simple, 5x5
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
    assert( mesh.coarse_boundary_cell(Point2(0.0, 2.0), 1) == 12 );
    assert( mesh.coarse_boundary_cell(Point2(0.0, 3.0), 4) == 12 );
    assert( mesh.coarse_boundary_cell(Point2(2.0, 5.0), 4) == 26 );
    assert( mesh.coarse_boundary_cell(Point2(4.0, 5.0), 3) == 27 );
    assert( mesh.coarse_boundary_cell(Point2(6.0, 4.0), 3) == 23 );
    assert( mesh.coarse_boundary_cell(Point2(6.0, 2.0), 2) == 17 );
    assert( mesh.coarse_boundary_cell(Point2(2.0, 0.0), 2) == 1 );
    assert( mesh.coarse_boundary_cell(Point2(3.5, 0.0), 1) == 3 );
    assert( mesh.coarse_boundary_cell(Point2(5.0, 0.0), 1) == 5 );


    // check surface indexing
    int s[2];
    Point2 p(0.5, 0.0);
    Point2 pcell(0.5, 0.5);
    int cell = mesh.coarse_cell_point( pcell );
    assert( cell == 0 );
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 65);

    // check a boundary corner on the south face
    p.x = 3.0; p.y = 0.0;
    pcell.x = 3.5; pcell.y = 0.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 3);
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 83 );
    // boundary corner on north face
    p.x = 4.0; p.y = 5.0;
    pcell.x = 4.5; pcell.y = 4.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 28);
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 94 );
    // boundary corner on east face
    p.x = 6.0; p.y = 4.0;
    pcell.x = 5.5; pcell.y = 3.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 23);
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 57 );
    // boundary corner on west face
    p.x = 0.0; p.y = 2.0;
    pcell.x = 0.5; pcell.y = 2.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 12);
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 44 );
    // interior x normal face
    p.x = 3.0; p.y = 1.5;
    pcell.x = 2.5; pcell.y = 1.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 8);
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 40 );
    // interior y normal face
    p.x = 2.5; p.y = 4.0;
    pcell.x = 2.5; pcell.y = 3.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 20);
    assert(mesh.coarse_surf_point( p, cell, s ) == 1 );
    assert(s[0] == 81 );
    // interior corner
    p.x = 4.0; p.y = 2.0;
    pcell.x = 4.5; pcell.y = 2.5;
    cell = mesh.coarse_cell_point( pcell );
    assert(cell == 16);
    assert(mesh.coarse_surf_point( p, cell, s ) == 2 );
    assert(s[0] == 48 );
    assert(s[1] == 85 );



    


    return 0;
}
