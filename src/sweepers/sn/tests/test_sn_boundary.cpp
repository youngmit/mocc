#include "UnitTest++/UnitTest++.h"

#include <algorithm>
#include <iostream>

#include "pugixml.hpp"

#include "core/angular_quadrature.hpp"
#include "core/constants.hpp"
#include "core/global_config.hpp"

#include "sn/sn_boundary.hpp"

using std::cout;
using std::endl;

using namespace mocc;

TEST( testboundary )
{
    // Start with the necessary boundary conditions, Mesh, and angular
    // quadrature to construct an Sn boundary condition.
    std::array<Boundary, 6> bc = {
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT
    };

    VecF x;
    for( real_t i=0.0; i<4.1; i++ ) {
        x.push_back(i);
    }
    VecF y;
    for( real_t i=0.0; i<5.1; i++ ) {
        y.push_back(i);
    }
    VecF z;
    for( real_t i=0.0; i<6.1; i++ ) {
        z.push_back(i);
    }

    // Should be a 4x5x6 mesh with unit cube regions
    Mesh mesh( 120, 120, x, y, z, bc );

    pugi::xml_document angquad_xml;
    pugi::xml_parse_result result =
        angquad_xml.load_string( "<ang_quad type=\"ls\" order=\"4\" />" );

    if( !result ) {
        cout << "failed to read ang quad" << endl;
    }

    AngularQuadrature ang_quad( angquad_xml.child("ang_quad") );

    // Actually make a BC object
    SnBoundary boundary( 1, ang_quad, mesh );
    SnBoundary out( 1, ang_quad, mesh );

    // Start testing things
    boundary.initialize( 3.14 );

    {
        ArrayF face = boundary.get_face( 0, 0, Normal::X_NORM );
        // Check face dimensions
        CHECK_EQUAL(face.size(), 5*6);
        // Check that the face is storing the right value
        for( int i=0; i<5*6; i++ ) {
            CHECK_EQUAL(3.14, face[i]);
        }
    }
    {
        ArrayF face = boundary.get_face( 0, 0, Normal::Y_NORM );
        CHECK_EQUAL(4*6, face.size());

        // Change the value of the face, set it, get it back and make sure it
        // stuck.
        face = boundary.get_face( 1, 1, Normal::X_NORM );
        face = 4.0;
        boundary.set_face( 1, 1, Normal::X_NORM, face );
        face = boundary.get_face( 1, 1, Normal::X_NORM );
        CHECK_EQUAL(4.0, face[3]);

        // Make sure that the rest of the faces are unaffected
        face = boundary.get_face( 0, 1, Normal::X_NORM );
        CHECK_EQUAL(3.14, face[3]);
    }
    {
        ArrayF face = boundary.get_face( 0, 0, Normal::Z_NORM );
        CHECK_EQUAL(4*5, face.size());
    }

    // Check the update routines. Stash some values in the out boundary, and use
    // them to update the boundary and make sure that the right values come back
    // out.
    {
        for( auto iface: AllNormals ) {
            for( int iang=0; iang<ang_quad.ndir(); iang++ ) {
                boundary.zero_face(0, iang, iface);
            }
        }
        ArrayF face = boundary.get_face( 0, 0, Normal::X_NORM );
        face = 1.77;
        out.set_face( 0, 0, Normal::X_NORM, face );

        int iangrefl = ang_quad.reflect(0, Normal::X_NORM);

        boundary.update( 0, 0, out );
        face = boundary.get_face( 0, iangrefl, Normal::X_NORM );
        for( int i=0; i<(int)face.size(); i++ ) {
            CHECK_EQUAL( 1.77, face[i] );
        }
    }
}

int main() {
    return UnitTest::RunAllTests();
}
