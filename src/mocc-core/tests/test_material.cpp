#include <iostream>
#include <string>
#include <cassert>

#include "global_config.hpp"
#include "file_scrubber.hpp"
#include "material_lib.hpp"
#include "material.hpp"
#include "fp_utils.hpp"

using std::cout;
using std::endl;

int main() {
    mocc::FileScrubber c5g7_file ("c5g7.xsl", "!");
    mocc::MaterialLib matlib(c5g7_file);

    matlib.assignID( 1, "MOX-4.3");

    const mocc::Material& mat = matlib.get_material_by_id(1);

    assert(mat.xsab().size() == 7);


    mocc::VecF out_scat { 1.702972340405E-01,
                          3.270915015982E-01,
                          4.558022000000E-01,
                          4.627124000000E-01,
                          2.862871691656E-01,
                          2.698171000000E-01,
                          2.735018000000E-01 };



    cout << endl << endl;
    for( int ig=0; ig<7; ig++ ) {
        assert(mocc::fp_equiv_ulp( mat.xssc().out(ig), out_scat[ig]) );


        
        const mocc::ScatRow& scat_row = mat.xssc().to(ig);
        unsigned int min_g = scat_row.min_g;
        unsigned int max_g = scat_row.max_g;
        for( int igg=min_g; igg<=max_g; igg++ ) {
            cout << scat_row.from[igg-min_g] << " ";
        }

        cout << endl;
    }
    assert(mat.xssc().to(3).from[0] == 5.04050E-09);
}
