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

    // Test the input/output
    bool isValidOutput() {
        H5Node h5file( "test_angquad.h5", H5Access::WRITE );
        ang_quad.output(h5file);
        AngularQuadrature new_ang_quad( h5file );
        return (new_ang_quad == ang_quad);
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

   
TEST_FIXTURE( LevelSymmetric_4, general )
{
    cout << ang_quad << endl;

    // Test the angle reflection capabilities
    test_reflect();

    // Test the angle reversal capabilities
    CHECK_EQUAL(7, ang_quad.reverse(1));
    CHECK_EQUAL(5, ang_quad.reverse(11));
    
    // Other tests
    CHECK_EQUAL(3, ang_quad.ndir_oct());
    // Test the weight sum is 8.0
    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
    // Test input/output
    CHECK(isValidOutput());
}

TEST_FIXTURE( LevelSymmetric_6, higher_order ) {
    
    // Test the weight sum is 8.0
    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
    // Test input/output
    CHECK(isValidOutput());
    
}


TEST_FIXTURE( AngQuadFixture, Chebyshev16Gauss3 ) {
    std::string inp = "<ang_quad type=\"cg\" azimuthal-order=\"16\" polar-order=\"3\" />";
    make_angquad( inp );
    // Test the weight sum is 8.0
   // CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);

}

/*
TEST_FIXTURE( AngQuadFixture, Chebyshev16Yamamoto3 ) {
    std::string inp = "<ang_quad type=\"cy\" azimuthal-order=\"16\" polar-order=\"3\" />";
    make_angquad( inp );
    // Test the weight sum is 8.0
//    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
}
*/
 
int main() {
    return UnitTest::RunAllTests();
}
