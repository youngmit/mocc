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

int main() {
    // For now, mostly just testing the coarse ray data, since its the easiers
    // to screw up.
    
    // Make a simple mesh, 1.0 cm pitch to keep things simple, 5x5
    mocc::VecF x;
    for( float xi=0.0; xi<5.1; xi++ ) {
        x.push_back(xi);
    }

    mocc::Mesh mesh( 25, 25, 5, 5, 1, x, x);

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
    
        for( int i=0; i< p.size(); i++ ) {
            cout << p[i] << endl;
    
            int cell = mesh.coarse_cell(cellpos[i]);
            cout << cell << endl;
            
            int s[2];
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
