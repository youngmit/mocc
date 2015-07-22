#include <string>
#include <iostream>
#include <cassert>

#include "pugixml.hpp"

#include "angular_quadrature.hpp"

using std::cout;
using std::endl;


int main() {
    std::string inp = "<ang_quad type=\"ls\" order=\"4\" />";

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string( inp.c_str() );

    if (!result) {
        cout << "failed to parse xml" << endl;   
    }


    mocc::AngularQuadrature ang_quad( doc.child("ang_quad") );

    assert(ang_quad.ndir_oct() == 3);

    // Test the angle reflection capabilities
    // octant 1
    assert(ang_quad.reflect(1, NORTH) == 10);
    assert(ang_quad.reflect(2, SOUTH) == 11);
    assert(ang_quad.reflect(2, EAST)  == 5);
    assert(ang_quad.reflect(0, WEST)  == 3);

    // octant 2
    assert(ang_quad.reflect(4, WEST)  == 1);
    assert(ang_quad.reflect(5, NORTH) == 8);

    // octant 3
    assert(ang_quad.reflect(7, WEST)  == 10);
    assert(ang_quad.reflect(6, SOUTH) == 3);

    // octant 4
    assert(ang_quad.reflect(11, EAST)  == 8);
    assert(ang_quad.reflect( 9, SOUTH) == 0);


    // Test the angle reversal capabilities
    assert(ang_quad.reverse(1) == 7);
    assert(ang_quad.reverse(11) == 5);
    

    return 0;
}