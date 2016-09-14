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

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "util/fp_utils.hpp"
#include "core/exponential.hpp"

using namespace mocc;

TEST(exp)
{

    Exponential_Linear<10000> exp;
    std::cout << "Max error: " << exp.max_error() << std::endl;

    real_t max_err = 0.0;
    for (real_t x = -10.0; x < 0.0; x += 0.1) {
        real_t exp_t = exp.exp(x);
        real_t exp_r = std::exp(x);
        max_err      = std::max(max_err, std::abs(exp_r - exp_t));
        std::cout << x << " " << exp_r << " " << exp_t << std::endl;
        CHECK(std::abs(exp_r - exp_t) < 2e-8);
    }

    std::cout << "max_err: " << max_err << std::endl;
}

TEST(exp_positive)
{
    Exponential_Linear<50000> exp(0.0, 10.0);

    std::cout << "max error from exp: " << exp.max_error() << std::endl;

    real_t max_err = 0.0;
    for (real_t x = 0.0; x < 10.0; x += 0.1) {
        real_t exp_t = exp.exp(x);
        real_t exp_r = std::exp(x);
        real_t err   = std::abs(exp_r - exp_t) / exp_r;
        max_err      = std::max(max_err, err);
        std::cout << x << " " << exp_r << " " << exp_t << " " << err
                  << std::endl;
        CHECK(err < 2e-8);
    }

    std::cout << "max_err: " << max_err << std::endl;
}

TEST(exp_coarse)
{
    Exponential_Linear<5> exp(-5.3, 0.0);

    CHECK_EQUAL(1.06, exp.dx());

    // Make sure the data points are right
    CHECK_CLOSE(std::exp(-5.3), exp[0], REAL_FUZZ);
    CHECK_CLOSE(std::exp(-4.24), exp[1], REAL_FUZZ);
    CHECK_CLOSE(std::exp(-3.18), exp[2], REAL_FUZZ);
    CHECK_CLOSE(std::exp(-2.12), exp[3], REAL_FUZZ);
    CHECK_CLOSE(std::exp(-1.06), exp[4], REAL_FUZZ);
    CHECK_CLOSE(std::exp(0.0), exp[5], REAL_FUZZ);

    // Spot-check some actual values
    // Start with the actual datapoints in the table
    CHECK_CLOSE(exp[0], exp.exp(-5.3) , REAL_FUZZ);
    CHECK_CLOSE(exp[1], exp.exp(-4.24) , REAL_FUZZ);
    CHECK_CLOSE(exp[2], exp.exp(-3.18) , REAL_FUZZ);
    CHECK_CLOSE(exp[3], exp.exp(-2.12) , REAL_FUZZ);
    CHECK_CLOSE(exp[4], exp.exp(-1.06) , REAL_FUZZ);
    CHECK_CLOSE(exp[5], exp.exp(0.0) , REAL_FUZZ);

    // Now some values away from the points
    CHECK_CLOSE(6.87479349415065E-03, exp.exp(-5.088), REAL_FUZZ);
    CHECK_CLOSE(7.29640444772866E-02, exp.exp(-2.756), REAL_FUZZ);

}

int main(int, const char *[])
{
    return UnitTest::RunAllTests();
}
