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
#include "fp_utils.hpp"

using namespace mocc;

TEST(fp_utils)
{
    double v1 = -1.249999999999996;
    double v2 = -1.25;
    std::cout << *((long*)&v1) - *((long*)&v2) << std::endl;
    CHECK(fp_equiv_saferel(v1, v2));
    CHECK(fp_equiv_saferel(-v1, -v2));
    CHECK(fp_equiv_saferel(v2, v1));
    CHECK(fp_equiv_saferel(-v2, -v1));

    v1 = 1.0;

    long i2 = *((long*)&v1) + 19;
    v2 = *((double*)&i2);
    std::cout << v1 << " " << v2 << " " << v2 - v1 << " " << std::numeric_limits<double>::epsilon() << std::endl;


}

int main()
{
    return UnitTest::RunAllTests();
}
