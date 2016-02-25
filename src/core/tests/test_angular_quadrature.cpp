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

class AngQuadFixture {
public: 
    AngularQuadrature ang_quad;
    
    void make_angquad( const std::string &inp_quad ) {
        pugi::xml_document doc;
        pugi::xml_parse_result result =  doc.load_string( inp_quad.c_str() );

        if (!result) {
            cout << "failed to parse xml" << endl;
        }

        ang_quad = AngularQuadrature( doc.child("ang_quad") );
    }

    void test_reflect() {
        int iang = 0;
        for( auto &angle: ang_quad ) {
            {
                int refl = ang_quad.reflect(iang, Normal::Y_NORM);
                auto &refl_ang = ang_quad[refl];

                CHECK_CLOSE(angle.ox, refl_ang.ox, 0.0000000000001);
                CHECK_CLOSE(angle.oy, -refl_ang.oy, 0.0000000000001);
                CHECK_CLOSE(angle.oz, refl_ang.oz, 0.0000000000001);
            }
            {
                int refl = ang_quad.reflect(iang, Normal::X_NORM);
                auto &refl_ang = ang_quad[refl];

                CHECK_CLOSE(angle.ox, -refl_ang.ox, 0.0000000000001);
                CHECK_CLOSE(angle.oy, refl_ang.oy, 0.0000000000001);
                CHECK_CLOSE(angle.oz, refl_ang.oz, 0.0000000000001);
            }
            {
                int refl = ang_quad.reflect(iang, Normal::Z_NORM);
                auto &refl_ang = ang_quad[refl];

                CHECK_CLOSE(angle.ox, refl_ang.ox, 0.0000000000001);
                CHECK_CLOSE(angle.oy, refl_ang.oy, 0.0000000000001);
                CHECK_CLOSE(angle.oz, -refl_ang.oz, 0.0000000000001);
            }

            iang++;
        }
    }

    real_t total_weight() {
        real_t wsum = 0.0;
        for( auto a: ang_quad ) {
            wsum += a.weight;
        }

        return wsum;
    }
};


class LevelSymmetric_4 : public AngQuadFixture {
public:
    LevelSymmetric_4() {
        this->make_angquad("<ang_quad type=\"ls\" order=\"4\" />");
    }
};


class LevelSymmetric_6 : public AngQuadFixture {
public:
    LevelSymmetric_6() {
        this->make_angquad("<ang_quad type=\"ls\" order=\"6\" />");
    }
};


/*

TEST_FIXTURE( AngQuadFixture, LevelSymmetric4 )
{
    std::string inp = "<ang_quad type=\"ls\" order=\"4\" />";
    make_angquad( inp );
    
    CHECK_EQUAL(3, ang_quad.ndir_oct());
    // Test the angle reflection capabilities
    // CHECK_EQUAL( expected, actual );
    
    // octant 1
    CHECK_EQUAL(10, ang_quad.reflect(1, Surface::NORTH));
    CHECK_EQUAL(11, ang_quad.reflect(2, Surface::SOUTH));
    CHECK_EQUAL(5, ang_quad.reflect(2, Surface::EAST));
    CHECK_EQUAL(3, ang_quad.reflect(0, Surface::WEST));

    // octant 2
    CHECK_EQUAL(1, ang_quad.reflect(4, Surface::WEST));
    CHECK_EQUAL(8, ang_quad.reflect(5, Surface::NORTH));

    // octant 3
    CHECK_EQUAL(10, ang_quad.reflect(7, Surface::WEST));
    CHECK_EQUAL(3, ang_quad.reflect(6, Surface::SOUTH));

    // octant 4
    CHECK_EQUAL(8, ang_quad.reflect(11, Surface::EAST));
    CHECK_EQUAL(0, ang_quad.reflect( 9, Surface::SOUTH));

    // octant 5
    CHECK_EQUAL(15, ang_quad.reflect(12, Surface::EAST));
    CHECK_EQUAL(17, ang_quad.reflect(14, Surface::EAST));
    CHECK_EQUAL(22, ang_quad.reflect(13, Surface::NORTH));



    // Test the angle reversal capabilities
    CHECK_EQUAL(7, ang_quad.reverse(1));
    CHECK_EQUAL(5, ang_quad.reverse(11));
}

TEST_FIXTURE( AngQuadFixture, LevelSymmetric6 ) {
    std::string inp = "<ang_quad type=\"ls\" order=\"6\" />";
    make_angquad( inp );
   
    // Test the weight sum is 8.0
    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
    
    // Test the input/output
    H5Node h5file("test_angquad.h5", H5Access::WRITE );
    ang_quad.output(h5file);
    AngularQuadrature new_ang_quad( h5file );
    CHECK(new_ang_quad == ang_quad);
}

*/

/*
TEST_FIXTURE( AngQuadFixture, Chebyshev16Gauss3 ) {
    std::string inp = "<ang_quad type=\"cg\" azimuthal-order=\"16\" polar-order=\"3\" />";
    make_angquad( inp );
    // Test the weight sum is 8.0
   // CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);

}
*/
/*
TEST_FIXTURE( AngQuadFixture, Chebyshev16Yamamoto3 ) {
    std::string inp = "<ang_quad type=\"cy\" azimuthal-order=\"16\" polar-order=\"3\" />";
    make_angquad( inp );
    // Test the weight sum is 8.0
//    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
}
*/
    
TEST_FIXTURE( LevelSymmetric_4, general )
{
    cout << ang_quad << endl;


    test_reflect();

    CHECK_EQUAL(3, ang_quad.ndir_oct());

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

/*    
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
*/
int main() {
    return UnitTest::RunAllTests();
}
