/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "UnitTest++/UnitTest++.h"

#include <iostream>
#include <vector>

#include "scattering_matrix.hpp"
#include "global_config.hpp"
#include "blitz_typedefs.hpp"

using std::cout;
using std::endl;

using namespace mocc;

const std::vector<VecF> sc  
    {
        {0.3,0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.3,0.0,0.0,0.0},
        {0.1,0.2,0.3,0.4,0.5,0.0,0.0},
        {0.0,0.1,0.2,0.3,0.4,0.5,0.0},
        {0.0,0.0,0.1,0.2,0.3,0.4,0.5},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.3}
    };

const std::vector<real_t> sc_dense 
    {
        0.3,0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.3,0.0,0.0,0.0,
        0.1,0.2,0.3,0.4,0.5,0.0,0.0,
        0.0,0.1,0.2,0.3,0.4,0.5,0.0,
        0.0,0.0,0.1,0.2,0.3,0.4,0.5,
        0.0,0.0,0.0,0.0,0.0,0.0,0.3
    };

// testing 1 group purely absorbing case
const std::vector<VecF> sc3
    { {0.0} };

ScatteringMatrix scat_matrix_ref(sc);

ArrayB2 vec_to_array( const std::vector<VecF> &vec ) {
    int ng = vec.size();
    ArrayB2 array(vec.size(),vec[0].size());
    for( int i=0; i<ng; i++ ) {
        for( int j=0; j<ng; j++ ) {
            array(i,j)=vec[i][j];
        }
    }
    return array;
};

//void printarray (real_t arg[], int length) {
//      for (int n=0; n<length; ++n)
//              cout << arg[n] << ' ';
//        cout << '\n';
//}

TEST( scat_matrix_vecorVecF ) {
    ScatteringMatrix scat_matrix(sc);
    int ng = sc.size();

    // test const ScatteringRow & to( int ig) const function
    ScatteringRow scat_row(0, 4, &sc[3][0]);
    CHECK( scat_row == scat_matrix.to(3) );

    // test copy constructor
    ScatteringMatrix scat_matrix_copy(scat_matrix);
    CHECK( scat_matrix_copy == scat_matrix );

    // test assignment operator
    ScatteringMatrix scat_matrix_assigned;
    scat_matrix_assigned = scat_matrix;
    CHECK( scat_matrix_assigned == scat_matrix );
    
    // test real_t self_scat (int group) const
    real_t self_scat_ref[ng];
    int ig = 0;
    for( auto &row : sc ) {
        self_scat_ref[ig] = row[ig];
        ig++;
    }
    
    real_t self_scat[ng];
    for( ig=0; ig<ng; ig++ ) {
        self_scat[ig] = scat_matrix.self_scat(ig);
    }

    //printarray(self_scat_ref,ng);
    //printarray(self_scat, ng);
    //CHECK_ARRAY_CLOSE( self_scat_ref, self_scat, ng, 0.0000000000001 );
    for( ig=0; ig<ng; ig++ ) {
        CHECK_CLOSE(self_scat_ref[ig], self_scat[ig], 0.0000000000001 );
    }
    //real_t bar [5] = { 10, 20, 30 };
    //real_t bar2 [5] = { 10, 20, 30 };
    //CHECK_ARRAY_CLOSE( bar, bar2, 5, 0.0000000000001 );
    
    // test n_group() function 
    CHECK_EQUAL( ng, scat_matrix.n_group() );
    
    // test out(unsigned int ig) const function
    real_t out_scat_ref[ng]={};
    
    for( auto &i : sc) {
        ig = 0;
        for( auto &j : i ) {
            out_scat_ref[ig] += j;
            ig++;
        }
    }

    for( ig=0; ig<ng; ig++ ) {
        CHECK_CLOSE(out_scat_ref[ig], scat_matrix.out(ig), 0.0000000000001 );
    }

    // test begin() and end() function
    ScatteringRow scat_row_2(0, 0, &sc[0][0]);
    CHECK( scat_row_2 == *(scat_matrix.begin()) );

    // test const ScatteringRow & to( int ig) const function
    ScatteringRow scat_row_3(6, 6, &sc[6][6]);
    CHECK( scat_row_3 == *(scat_matrix.end()-1) );
    // followup test of begin() and end() methods of ScatteringRow
    CHECK_CLOSE( 0.3, *(scat_row_3.begin()), 0.0000000000001 );
    CHECK_CLOSE( 0.3, *(scat_row_3.end()-1), 0.0000000000001 );
    
    //test VecF as_vector() member function
    CHECK( sc_dense == scat_matrix.as_vector() ); 
    
    // test operator== 
    CHECK( scat_matrix_ref == scat_matrix );
    
    // test operator!=
    std::vector<VecF> sc2(sc);
    sc2[1][1]=0.3;
    ScatteringMatrix scat_matrix2(sc2);
    CHECK( scat_matrix != scat_matrix2 );

    // test indexing into a specific element scat(to, from)
    CHECK_CLOSE( 0.0, scat_matrix.to(1)[1], 0.0000000000001 );
    CHECK_CLOSE( 0.4, scat_matrix.to(3)[3], 0.0000000000001 );

}

TEST( scat_matrix_ArrayB2 ) {
    ScatteringMatrix scat_matrix(vec_to_array(sc));
    // As long as TEST scat_matrix_vectorVecF passed, showing scat_matrix in
    // both tests are equal to scat_matrix_ref will surfice. 
    CHECK(scat_matrix_ref==scat_matrix);

}

TEST( vecF_purely_absorbing ) {
    ScatteringMatrix scat_matrix(sc3);
    
    CHECK_CLOSE( 0.0, scat_matrix.to(0)[0], 0.0000000000001 );
}

int main() {
    return UnitTest::RunAllTests();
}
