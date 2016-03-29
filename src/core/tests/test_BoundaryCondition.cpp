#include "UnitTest++/UnitTest++.h"

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "boundary_condition.hpp"

using namespace mocc;
using std::cout;
using std::endl;

// This is a fixture that provides an in and out boundary condition that look a
// lot like you would expect for an MoC sweeper. there are only boundary values
// on the X and Y normal faces, and a different number on each, depending on
// angle
class BCIrregularFixture {
    public:
        BCIrregularFixture():
            angquad(make_angquad()),
            nang(angquad.ndir()/2),
            ngroup(2),
            nbc(
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
                    }),
            bc_(
                {
                    Boundary::REFLECT,
                    Boundary::REFLECT,
                    Boundary::REFLECT,
                    Boundary::REFLECT,
                    Boundary::REFLECT,
                    Boundary::REFLECT
                }),
            in(2, angquad, bc_, nbc),
            out(1, angquad, bc_, nbc)
        {
            return;
        }

    // Data
    public:
        AngularQuadrature angquad;
        int nang;
        int ngroup;
        std::vector< BC_Size_t > nbc;
    private:
        BC_Type_t bc_;
    public:
        BoundaryCondition in;
        BoundaryCondition out;

    private:
        AngularQuadrature make_angquad() {
            std::string angquad_xml = "<ang_quad type=\"ls\" order=\"4\" />";
            pugi::xml_document angquad_node;
            angquad_node.load_string(angquad_xml.c_str());
            return AngularQuadrature(angquad_node.child("ang_quad"));
        }
};

TEST_FIXTURE( BCIrregularFixture, test_bc )
{
    CHECK_EQUAL( 192, in.size() );
    
    in.initialize_scalar( 7.345 );
    for( int ig=0; ig<ngroup; ig++ ) {
        for( int ia=0; ia<nang; ia++ ) {
            auto face = in.get_face( ig, ia, Normal::X_NORM );
            CHECK_EQUAL(nbc[ia][(int)Normal::X_NORM], face.first);
            for( int ibc=0; ibc<8; ibc++ ) {
                CHECK_EQUAL( 7.345, face.second[ibc] );
            }
        }
    }

    ArrayB1 spectrum(2);
    spectrum(0) = 2.2222;
    spectrum(1) = 4.4444;
    in.initialize_spectrum( spectrum );
    for( int ia=0; ia<nang; ia++ ) {
        auto face = in.get_face( 0, ia, Normal::X_NORM );
        CHECK_EQUAL(nbc[ia][(int)Normal::X_NORM], face.first);
        for( int ibc=0; ibc<8; ibc++ ) {
            CHECK_EQUAL( 2.2222, face.second[ibc] );
        }
    }
    for( int ia=0; ia<nang; ia++ ) {
        auto face = in.get_face( 1, ia, Normal::X_NORM );
        CHECK_EQUAL(nbc[ia][(int)Normal::X_NORM], face.first);
        for( int ibc=0; ibc<8; ibc++ ) {
            CHECK_EQUAL( 4.4444, face.second[ibc] );
        }
    }

    out.initialize_scalar(0.0);

    auto outface = out.get_face(0, 0, Normal::X_NORM);
    for( int i=0; i<nbc[0][(int)Normal::X_NORM]; i++ ) {
        outface.second[i] = 3.3333;
    }

    in.update(0, out);

    int iang_refl = angquad.reflect(0, Normal::X_NORM);
    auto inface = in.get_face(0, iang_refl, Normal::X_NORM);
    for( int ibc=0; ibc<5; ibc++ ) {
        CHECK_EQUAL(3.3333, inface.second[ibc]);
    }
    for( int ibc=5; ibc<8; ibc++ ) {
        CHECK_EQUAL(0.0, inface.second[ibc]);
    }

    std::cout << in << std::endl;

}

int main() {
    return UnitTest::RunAllTests();
}
