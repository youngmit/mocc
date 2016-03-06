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

class ChebyshevGauss_16_3 : public AngQuadFixture {
public:
    ChebyshevGauss_16_3() {
        this->make_angquad("<ang_quad type=\"cg\" azimuthal-order=\"16\" polar-order=\"3\" />");
    }
};

class ChebyshevYamamoto_16_3 : public AngQuadFixture {
public:
    ChebyshevYamamoto_16_3() {
        this->make_angquad("<ang_quad type=\"cy\" azimuthal-order=\"16\" polar-order=\"3\" />");
    }
};
   
TEST_FIXTURE( LevelSymmetric_4, general )
{
    //cout << ang_quad << endl;

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
    
    // Test the angle reflection capabilities
    test_reflect();

    // Other tests
    CHECK_EQUAL(6, ang_quad.ndir_oct());
    // Test the weight sum is 8.0
    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
    // Test input/output
    CHECK(isValidOutput());
}


TEST_FIXTURE( ChebyshevYamamoto_16_3, cy_general ) {
    std::string inp = "<ang_quad type=\"cy\" azimuthal-order=\"16\" polar-order=\"3\" />";
    // Test the angle reflection capabilities
    test_reflect();

    // Other tests
    CHECK_EQUAL(48, ang_quad.ndir_oct());
    // Test the weight sum is 8.0
    CHECK_CLOSE(8.0, total_weight(), 0.00000000000001);
/**
 * \todo isValidOutput() fails when Angle constructor used to generate the
 * output angular quadrature is different from the Angle constructor used in the
 * constructor that accepts a h5file node.
 */
    // Test input/output
    // CHECK(isValidOutput());
    // Test the first angle
    CHECK_CLOSE(0.049087385212340, ang_quad[0].alpha,    0.0000000000001); 
    CHECK_CLOSE(0.167429147795000, ang_quad[0].theta,    0.0000000000001); 
    CHECK_CLOSE(0.166447265186000, ang_quad[0].ox,       0.0000000000001); 
    CHECK_CLOSE(0.008177029791330, ang_quad[0].oy,       0.0000000000001); 
    CHECK_CLOSE(0.986016452244020, ang_quad[0].oz,       0.0000000000001); 
    CHECK_CLOSE(0.002889562500000, ang_quad[0].weight,   0.0000000000001); 
    CHECK_CLOSE(6.000672075260800, ang_quad[0].rsintheta,0.0000000000001); 
}


TEST_FIXTURE( ChebyshevGauss_16_3, cg_general ) {
    
    //cout << ang_quad << endl;
    // Test the angle reflection capabilities
    test_reflect();

    // Other tests
    CHECK_EQUAL(48, ang_quad.ndir_oct());
    // Test the weight sum is 8.0
    CHECK_CLOSE(8.0, total_weight(), 0.0000000000001);
    // Test the first angle
    CHECK_CLOSE(0.049087385212340, ang_quad[0].alpha,    0.000000000001); 
    CHECK_CLOSE(0.374822141002272, ang_quad[0].theta,    0.000000000001); 
    CHECK_CLOSE(0.365666032196437, ang_quad[0].ox,       0.000000000001); 
    CHECK_CLOSE(0.017964020229512, ang_quad[0].oy,       0.000000000001); 
    CHECK_CLOSE(0.930572752059133, ang_quad[0].oz,       0.000000000001); 
    CHECK_CLOSE(0.029244620910793, ang_quad[0].weight,   0.000000000001); 
    CHECK_CLOSE(2.731441720757410, ang_quad[0].rsintheta,0.000000000001); 
  
    // Test input/output
    // CHECK(isValidOutput());
}


int main() {
    return UnitTest::RunAllTests();
}






