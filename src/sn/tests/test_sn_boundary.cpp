#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <iostream>

#include "constants.hpp"
#include "global_config.hpp"
#include "sn_boundary.hpp"

using std::cout;
using std::endl;

using namespace mocc;

BOOST_AUTO_TEST_CASE( testboundary )
{
    // SnBoundary( int n_grp, int n_ang , int nx, int ny, int nz);
    SnBoundary boundary( 2, 8, 4, 5, 6 );

    boundary.initialize( 3.14 );

    {
        ArrayF face = boundary.get_face( 0, 0, Normal::X_NORM );
        BOOST_CHECK(face.size() == 5*6);
        BOOST_CHECK(face[0] == 3.14);
    }
    {
        ArrayF face = boundary.get_face( 0, 0, Normal::Y_NORM );
        BOOST_CHECK(face.size() == 4*6);

        // Change the value of the face, set it, get it back and make sure it
        // stuck.
        face = boundary.get_face( 1, 1, Normal::X_NORM );
        face = 4.0;
        boundary.set_face( 1, 1, Normal::X_NORM, face );
        face = boundary.get_face( 1, 1, Normal::X_NORM );
        BOOST_CHECK(face[3] == 4.0);

        // Make sure that the rest of the faces are unaffected
        face = boundary.get_face( 0, 1, Normal::X_NORM );
        BOOST_CHECK(face[3] == 3.14);
    }
    {
        ArrayF face = boundary.get_face( 0, 0, Normal::Z_NORM );
        BOOST_CHECK(face.size() == 4*5);
    }
}
