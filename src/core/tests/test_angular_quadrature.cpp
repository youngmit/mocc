#undef __STRICT_ANSI__
#undef _REENT_ONLY
#include <stdlib.h>
#include <string>
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <cassert>

#include "pugixml.hpp"

#include "angular_quadrature.hpp"

using std::cout;
using std::endl;

using namespace mocc;

BOOST_AUTO_TEST_CASE( testall )
{
    std::string inp = "<ang_quad type=\"ls\" order=\"4\" />";

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string( inp.c_str() );

    if (!result) {
        cout << "failed to parse xml" << endl;
    }


    mocc::AngularQuadrature ang_quad( doc.child("ang_quad") );

    cout << ang_quad << endl;

    BOOST_CHECK(ang_quad.ndir_oct() == 3);

    // Test the angle reflection capabilities
    // octant 1
    BOOST_CHECK_EQUAL(ang_quad.reflect(1, Surface::NORTH), 10);
    BOOST_CHECK_EQUAL(ang_quad.reflect(2, Surface::SOUTH), 11);
    BOOST_CHECK_EQUAL(ang_quad.reflect(2, Surface::EAST), 5);
    BOOST_CHECK_EQUAL(ang_quad.reflect(0, Surface::WEST), 3);

    // octant 2
    BOOST_CHECK_EQUAL(ang_quad.reflect(4, Surface::WEST), 1);
    BOOST_CHECK_EQUAL(ang_quad.reflect(5, Surface::NORTH), 8);

    // octant 3
    BOOST_CHECK_EQUAL(ang_quad.reflect(7, Surface::WEST), 10);
    BOOST_CHECK_EQUAL(ang_quad.reflect(6, Surface::SOUTH), 3);

    // octant 4
    BOOST_CHECK_EQUAL(ang_quad.reflect(11, Surface::EAST), 8);
    BOOST_CHECK_EQUAL(ang_quad.reflect( 9, Surface::SOUTH), 0);

    // octant 5
    BOOST_CHECK_EQUAL(ang_quad.reflect(12, Surface::EAST), 15);
    BOOST_CHECK_EQUAL(ang_quad.reflect(14, Surface::EAST), 17);
    BOOST_CHECK_EQUAL(ang_quad.reflect(13, Surface::NORTH), 22);



    // Test the angle reversal capabilities
    BOOST_CHECK_EQUAL(ang_quad.reverse(1), 7);
    BOOST_CHECK_EQUAL(ang_quad.reverse(11), 5);

}
