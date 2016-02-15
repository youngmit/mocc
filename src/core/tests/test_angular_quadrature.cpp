#include "UnitTest++/UnitTest++.h"

#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "pugixml.hpp"

#include "angular_quadrature.hpp"

using std::cout;
using std::endl;

using namespace mocc;

class LevelSymmetric_4 {
private:
    static AngularQuadrature make_angquad() {
        std::string inp = "<ang_quad type=\"ls\" order=\"4\" />";
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_string( inp.c_str() );

        if (!result) {
            cout << "failed to parse xml" << endl;
        }

        return AngularQuadrature( doc.child("ang_quad") );
    }
public:
    AngularQuadrature ang_quad;
    LevelSymmetric_4(): ang_quad(LevelSymmetric_4::make_angquad()) {
        return;
    }
};

class LevelSymmetric_6 {
private:
    static AngularQuadrature make_angquad() {
        std::string inp = "<ang_quad type=\"ls\" order=\"6\" />";
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_string( inp.c_str() );

        if (!result) {
            cout << "failed to parse xml" << endl;
        }

        return AngularQuadrature( doc.child("ang_quad") );
    }
public:
    AngularQuadrature ang_quad;
    LevelSymmetric_6(): ang_quad(LevelSymmetric_6::make_angquad()) {
        return;
    }
};

TEST_FIXTURE( LevelSymmetric_4, general )
{
    cout << ang_quad << endl;

    CHECK_EQUAL(ang_quad.ndir_oct(), 3);

    // Test the angle reflection capabilities
    // octant 1
    CHECK_EQUAL(ang_quad.reflect(1, Surface::NORTH), 10);
    CHECK_EQUAL(ang_quad.reflect(2, Surface::SOUTH), 11);
    CHECK_EQUAL(ang_quad.reflect(2, Surface::EAST), 5);
    CHECK_EQUAL(ang_quad.reflect(0, Surface::WEST), 3);

    // octant 2
    CHECK_EQUAL(ang_quad.reflect(4, Surface::WEST), 1);
    CHECK_EQUAL(ang_quad.reflect(5, Surface::NORTH), 8);

    // octant 3
    CHECK_EQUAL(ang_quad.reflect(7, Surface::WEST), 10);
    CHECK_EQUAL(ang_quad.reflect(6, Surface::SOUTH), 3);

    // octant 4
    CHECK_EQUAL(ang_quad.reflect(11, Surface::EAST), 8);
    CHECK_EQUAL(ang_quad.reflect( 9, Surface::SOUTH), 0);

    // octant 5
    CHECK_EQUAL(ang_quad.reflect(12, Surface::EAST), 15);
    CHECK_EQUAL(ang_quad.reflect(14, Surface::EAST), 17);
    CHECK_EQUAL(ang_quad.reflect(13, Surface::NORTH), 22);



    // Test the angle reversal capabilities
    CHECK_EQUAL(ang_quad.reverse(1), 7);
    CHECK_EQUAL(ang_quad.reverse(11), 5);
}

TEST_FIXTURE( LevelSymmetric_6, higher_order ) {
    real_t wsum = 0.0;

    for( auto a: ang_quad ) {
        wsum += a.weight;
    }

    CHECK_CLOSE(8.0, wsum, 0.00000000000001);
}

TEST_FIXTURE( LevelSymmetric_6, input_output ) {
    H5Node h5file("test_angquad.h5", H5Access::WRITE );
    ang_quad.output(h5file);

    AngularQuadrature new_ang_quad( h5file );

    CHECK(new_ang_quad == ang_quad);
}

int main() {
    return UnitTest::RunAllTests();
}
