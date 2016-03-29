#include "UnitTest++/UnitTest++.h"

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
using namespace mocc;

TEST( material )
{
    FileScrubber c5g7_file("c5g7.xsl", "!");
    MaterialLib matlib(c5g7_file);

    matlib.assignID( 1, "MOX-4.3" );

    const Material& mat = matlib.get_material_by_id(1);

    CHECK_EQUAL(7, mat.xsab().size());

    VecF out_scat { 1.702972340405E-01,
                          3.270915015982E-01,
                          4.558022000000E-01,
                          4.627124000000E-01,
                          2.862871691656E-01,
                          2.698171000000E-01,
                          2.735018000000E-01 };
    VecI min_g { 0, 0, 0, 0, 3, 4, 4 };
    VecI max_g { 0, 1, 2, 4, 5, 6, 6 };
    



    for( int ig=0; ig<7; ig++ ) {
        CHECK_CLOSE( mat.xssc().out(ig), out_scat[ig], 0.000000000001 );

        const ScatteringRow& scat_row = mat.xssc().to(ig);
        CHECK_EQUAL( min_g[ig], scat_row.min_g );
        CHECK_EQUAL( max_g[ig], scat_row.max_g );
    }
    CHECK_CLOSE(5.04050E-09, mat.xssc().to(3).from[0], 0.000000000001);
}

int main(int, const char*[]) {
    return UnitTest::RunAllTests();
}
