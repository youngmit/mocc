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
    // For now, mostly just testing the coarse ray data, since its the easiers
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

    // Everything below should really be in a test for the ray types, rather
    // than the mesh test.

    // Make a nasty couple of points, should start on a corner of the mesh on
    // tha boundary, pass through some corners on the interior and exit on a
    // corner on the east boundary
    
    {
        mocc::Point2 p1( 3.0, 0.0 );
        mocc::Point2 p2( 5.0, 2.0 );
    
        std::vector<mocc::Point2> p;
        p.push_back(p1);
        p.push_back(p2);
    
        mesh.trace(p);
    
        mocc::Position cellpos[3] = { 
                            mocc::Position(3, 0, 0), 
                            mocc::Position(3, 0, 0), 
                            mocc::Position(4, 1, 0) };
        mocc::VecI surf;
    
        for( unsigned int i=0; i< p.size(); i++ ) {
            cout << p[i] << endl;
    
            cell = mesh.coarse_cell(cellpos[i]);
            cout << cell << endl;
            
            int nsurf = mesh.coarse_surf_point(p[i], cell, s);
            for( int i=0; i<nsurf; i++ ) {
                surf.push_back(s[i]);
            }
        }
    
        cout << "surface crossings:" << endl;
    
        for( auto &si: surf ) {
            cout << si << endl;
        }
        cout << endl;
    }

    // Try a simpler ray
    {
        mocc::Point2 p1( 1.2, 0.0 );
        mocc::Point2 p2( 2.8, 5.0 );
        std::vector<mocc::Point2> p;
        p.push_back(p1);
        p.push_back(p2);

        mesh.trace(p);

        mocc::VecI surf;
        int s[2];


        auto p_prev = *(p.begin());
        int cell = mesh.coarse_cell_point(mocc::Midpoint(p[0], p[1]));
        int nsurf = mesh.coarse_surf_point(p_prev, cell, s);
        for( int i=0; i<nsurf; i++ ) {
            surf.push_back(s[i]);
        }
        for( auto pi=p.begin()+1; pi!=p.end(); ++pi ) {
            cell = mesh.coarse_cell_point(mocc::Midpoint(*pi, p_prev));
            int nsurf = mesh.coarse_surf_point(*pi, cell, s);
            for( int i=0; i<nsurf; i++ ) {
                surf.push_back(s[i]);
            }
            p_prev = *pi;
        }

        for( auto &si: surf ) {
            cout << si << endl;
        }


    }

    return 0;
}
