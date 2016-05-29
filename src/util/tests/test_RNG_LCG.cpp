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

#include "rng_lcg.hpp"

#include <fstream>
#include <iostream>


TEST(test_RNG_LCG) {
    {
        mocc::RNG_LCG rng;
        std::ofstream fout("rand1.txt");

        for( int i=0; i<10000; i++ ) {
            fout << rng.random() << std::endl;
        }
    }
    {
        mocc::RNG_LCG rng(1);
        std::ofstream fout("rand2.txt");

        for( int i=0; i<10000; i++ ) {
            fout << rng.random() << std::endl;
        }
    }
}

TEST(uniformity) {
    mocc::RNG_LCG rng;
    std::vector<int> histogram(100, 0);
    int N = 100000000;
    for(int i=0; i<N; i++) {
        histogram[int(rng.random()*100)]++;
    }

    std::ofstream fout("histogram_uniform.txt");
    for( const auto &v: histogram ) {
        double vs = (double)v/(N/100);
        std::cout << vs << std::endl;
        fout << vs << std::endl;
    }
}

int main() {
    return UnitTest::RunAllTests();
}
