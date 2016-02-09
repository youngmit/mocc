#include "UnitTest++/UnitTest++.h"

#include <iostream>
#include <string>

#include "pugixml.hpp"

#include "core/blitz_typedefs.hpp"
#include "core/core_mesh.hpp"
#include "core/eigen_interface.hpp"
#include "core/global_config.hpp"
#include "core/material_lib.hpp"

#include "sweepers/moc/moc_sweeper.hpp"

using namespace mocc;
using std::cout;
using std::endl;
using std::cin;
using moc::MoCSweeper;

// This unit test makes sure that the MoC sweeper is capable of solving an
// infinite homogeneous medium case. The actual geometry isnt square, to try and
// catch certain errors. We do a little Eigen linear system solve to get the
// "right" answer to compare against.
//

std::string ihm_xml = 
"<mesh id=\"1\" type=\"rect\" pitch=\"1.26\">"
"<sub_x>3</sub_x>"
"<sub_y>3</sub_y>"
"</mesh>"
"<pin id=\"1\" mesh=\"1\">"
"1 1 1 1 1 1 1 1 1"
"</pin>"
"<lattice id=\"1\" nx=\"3\" ny=\"2\">"
"1 1 1 1 1 1"
"</lattice>"
"<assembly id=\"1\" np=\"1\" hz=\"0.5\">"
"<lattices>"
"1"
"</lattices>"
"</assembly>"
"<core nx=\"1\" ny=\"1\""
" north=\"reflect\""
" south=\"reflect\""
" east=\"reflect\""
" west=\"reflect\""
" top=\"reflect\""
" bottom=\"reflect\" >"
"1"
"</core>"
""
"<material_lib path=\"c5g7.xsl\">"
"<material id=\"1\" name=\"UO2-3.3\" />"
"</material_lib>"
""
"<source scattering=\"P0\" />"
"<sweeper type=\"moc\" n_inner=\"10\">"
"    <ang_quad type=\"ls\" order=\"2\" />"
"    <rays spacing=\"0.01\" />"
"</sweeper>";

pugi::xml_document xml_doc;

// Basic extension of the MoCSweeper class, which allows us to set a flux
// spectrum
class TestMoCSweeper : public MoCSweeper {
public:
    TestMoCSweeper( const pugi::xml_node &input, const CoreMesh &mesh ):
        MoCSweeper( input, mesh ) 
    {
        return;
    }

    void set_spectrum( const ArrayB1 &spectrum ) {
        int ig = 0;
        for( auto v: spectrum ) {
            flux_(blitz::Range::all(), ig) = v;
            ig++;
        }
        return;
    }
};

// This routine generates the reference solution
void reference_solution(real_t &k_eff, ArrayB1 &flux, ArrayB1 &psi );

TEST( moc_ihm )
{
    auto result = xml_doc.load_string(ihm_xml.c_str());
    CHECK(result);
    if(!result) {
        cout << result.description() << endl;
        cout << result.offset << endl;
        cout << "\"" << ihm_xml.substr(result.offset-10, 20) << "\"" << endl;
    }

    int ng = 7;
    ArrayB1 flux_ref(ng);
    ArrayB1 psi_ref(ng);
    real_t k_ref;
    reference_solution( k_ref, flux_ref, psi_ref );
    cout << "reference flux: " << flux_ref << endl;


    CoreMesh core_mesh(xml_doc);

    TestMoCSweeper sweeper(xml_doc.child("sweeper"), core_mesh);
    auto source = sweeper.create_source(xml_doc.child("source"));
    sweeper.assign_source(source.get());

    // Set the flux on the sweeper to match the true solution spectrum.
    sweeper.set_spectrum(flux_ref);

    // Calculate fission source for the true solution
    ArrayB1 fission_source(sweeper.n_reg());
    fission_source = 0.0;
    sweeper.calc_fission_source( k_ref, fission_source );

    for( int ireg=0; ireg<sweeper.n_reg(); ireg++ ) {
        CHECK_CLOSE(1.0, fission_source(ireg), 0.000000000000001);
    }

    for( int ig=0; ig<ng; ig++ ) {
        cout << "pre sweep: "<< sweeper.flux()(blitz::Range::all(),ig) << endl;
        source->initialize_group(ig);
        source->fission(fission_source, ig);
        source->in_scatter(ig);
        sweeper.sweep(ig);
        cout << sweeper.flux()(blitz::Range::all(),ig) << endl;
    }

}

void reference_solution( real_t &k_eff, ArrayB1 &flux, ArrayB1 &psi ) {
    const MaterialLib mat_lib(xml_doc.child("material_lib"));
    const Material &mat = mat_lib[1];

    Eigen::Matrix<real_t, 7, 7> M;
    Eigen::Matrix<real_t, 7, 1> chi;
    Eigen::Matrix<real_t, 1, 7> nf;

    M.fill(0.0);
    for( int ig=0; ig<mat_lib.n_group(); ig++ ) {
        M(ig, ig) += mat.xstr(ig);
        const ScatteringRow &scat_row = mat.xssc().to(ig);
        for( int igg=scat_row.min_g; igg<=scat_row.max_g; igg++ ) {
            M(ig, igg) -= scat_row[igg];
        }
        chi(ig,0) = mat.xsch(ig);
        nf(0,ig) = mat.xsnf(ig);
    }
    cout << "M: " << M << endl;

    Eigen::Matrix<real_t, 7, 1> phi = M.inverse()*chi;
    k_eff = nf*phi;

    real_t psi_tot = nf*phi;
    auto psi_eigen = psi_tot*chi/k_eff;

    for( int ig=0; ig<mat_lib.n_group(); ig++ ) {
        flux(ig) = phi(ig, 0);
        psi(ig) = psi_eigen(ig, 0);
    }

    return;
}

int main() {
    return UnitTest::RunAllTests();
}
