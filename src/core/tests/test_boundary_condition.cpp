#include "UnitTest++/UnitTest++.h"

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "boundary_condition.hpp"

using namespace mocc;
using std::cout;
using std::endl;

TEST( test_bc )
{
    std::string angquad_xml = "<ang_quad type=\"ls\" order=\"4\" />";
    pugi::xml_document angquad_node;
    angquad_node.load_string(angquad_xml.c_str());
    AngularQuadrature angquad(angquad_node.child("ang_quad"));

    int nang = 12;
    int ngroup = 2;
    std::vector< BC_Size_t > nbc = 
    {
        {5, 3, 0},
        {4, 4, 0},
        {3, 5, 0},
        {5, 3, 0},
        {4, 4, 0},
        {3, 5, 0},
        {5, 3, 0},
        {4, 4, 0},
        {3, 5, 0},
        {5, 3, 0},
        {4, 4, 0},
        {3, 5, 0}
    };

    BC_Type_t bc =
    {
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::REFLECT,
        Boundary::VACUUM,
        Boundary::REFLECT,
        Boundary::REFLECT
    };

    BoundaryCondition in( ngroup, angquad, bc, nbc );

    CHECK_EQUAL( 192, in.size() );
    
    in.initialize_scalar( 7.345 );
    for( int ig=0; ig<ngroup; ig++ ) {
        for( int ia=0; ia<nang; ia++ ) {
            const real_t *face = in.get_face( ig, ia, Normal::X_NORM );
            for( int ibc=0; ibc<8; ibc++ ) {
                CHECK_EQUAL( 7.345, face[ibc] );
            }
        }
    }

    ArrayB1 spectrum(2);
    spectrum(0) = 2.2222;
    spectrum(1) = 4.4444;
    in.initialize_spectrum( spectrum );
    for( int ia=0; ia<nang; ia++ ) {
        const real_t *face = in.get_face( 0, ia, Normal::X_NORM );
        for( int ibc=0; ibc<8; ibc++ ) {
            CHECK_EQUAL( 2.2222, face[ibc] );
        }
    }
    for( int ia=0; ia<nang; ia++ ) {
        const real_t *face = in.get_face( 1, ia, Normal::X_NORM );
        for( int ibc=0; ibc<8; ibc++ ) {
            CHECK_EQUAL( 4.4444, face[ibc] );
        }
    }

    BoundaryCondition out( 1, angquad, bc, nbc );
    out.initialize_scalar(0.0);

    real_t *outface = out.get_face(0, 0, Normal::X_NORM);
    for( int i=0; i<nbc[0][(int)Normal::X_NORM]; i++ ) {
        outface[i] = 3.3333;
    }

    in.update(0, out);

    int iang_refl = angquad.reflect(0, Normal::X_NORM);
    real_t *inface = in.get_face(0, iang_refl, Normal::X_NORM);
    for( int ibc=0; ibc<5; ibc++ ) {
        CHECK_EQUAL(3.3333, inface[ibc]);
    }
    for( int ibc=5; ibc<8; ibc++ ) {
        CHECK_EQUAL(0.0, inface[ibc]);
    }
}

int main(int, const char*[]) {
    return UnitTest::RunAllTests();
}
