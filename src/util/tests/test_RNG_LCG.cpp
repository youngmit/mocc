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

#include "fp_utils.hpp"
#include "rng_lcg.hpp"

#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

using namespace mocc;

TEST(test_RNG_LCG)
{
    {
        mocc::RNG_LCG rng;
        std::ofstream fout("rand1.txt");

        for (int i = 0; i < 10000; i++) {
            fout << rng.random() << std::endl;
        }
    }
    {
        mocc::RNG_LCG rng(1);
        std::ofstream fout("rand2.txt");

        for (int i = 0; i < 10000; i++) {
            fout << rng.random() << " " << rng.random() << std::endl;
        }
        real_t v = rng.random();
        cout << "current random value: " << v << endl;
        // reset and skip ahead
        rng.set_seed(1);
        rng.jump_ahead(20000);
        CHECK_EQUAL(v, rng.random());
    }
}

TEST(uniformity_standard)
{
    mocc::RNG_LCG rng;
    std::vector<int> histogram(100, 0);
    int N = 10000000;
    for (int i = 0; i < N; i++) {
        histogram[int(rng.random() * 100)]++;
    }

    std::ofstream fout("histogram_uniform.txt");
    real_t max_diff = 0.0;
    real_t rms_diff = 0.0;
    for (const auto &v : histogram) {
        double vs = (double)v / (N / 100);
        max_diff  = std::max(max_diff, std::abs(vs - 1.0));
        rms_diff += (vs - 1.0) * (vs - 1.0);
        fout << vs << std::endl;
    }
    std::cout << "Max variation from 1: " << max_diff << std::endl;
    CHECK(max_diff < 0.008);

    std::cout << "RMS: " << rms_diff << std::endl;
}

TEST(uniformity_custom_bounds)
{
    mocc::RNG_LCG rng;
    std::vector<int> histogram(100, 0);
    int N = 1000000;
    for (int i = 0; i < N; i++) {
        // cout << rng.random(-5.0, 1.0) << endl;
        // std::cin.ignore();
        //        histogram[int(rng.random() * 100)]++;
    }
}

TEST(estimate_pi)
{
    mocc::RNG_LCG rng;
    int n_in = 0;

    int N = 10000000;
    for (int i = 0; i < N; i++) {
        real_t x = rng.random();
        real_t y = rng.random();
        if (x * x + y * y < 1.0) {
            n_in++;
        }
    }
    std::cout << "pi estimate: " << 4.0 * n_in / N << std::endl;
}

TEST(random_bitstream)
{
    // this makes a 1MB binary file of pseudorandomness, which can be fed to
    // external tools to test it. Keep in mind, however, that the RNG only makes
    // 63 bits of randomness, so every 64th bit in the file is deterninistically
    // zero, which will result in poor results from some tests
    mocc::RNG_LCG rng;
    std::ofstream fout("random.binary", std::ios::binary | std::ios::out);
    std::vector<unsigned long> stream;
    int N = (1 << 20) / sizeof(unsigned long);
    stream.reserve(N);
    for (int i = 0; i < N; i++) {
        stream.push_back(rng());
    }

    fout.write((char *)stream.data(),
               (std::streamsize)N * sizeof(unsigned long));
}

TEST(sample_cdf)
{
    mocc::RNG_LCG rng;

    VecF pdf;
    pdf.push_back(0.05);
    pdf.push_back(0.075);
    pdf.push_back(0.1);
    pdf.push_back(0.125);
    pdf.push_back(0.1);
    pdf.push_back(0.125);
    pdf.push_back(0.15);
    pdf.push_back(0.075);
    pdf.push_back(0.1);
    pdf.push_back(0.1);

    VecF cdf;
    real_t prev = 0.0;
    for (const auto &p : pdf) {
        cdf.push_back(prev + p);
        prev = cdf.back();
    }
    CHECK_CLOSE(1.0, cdf.back(), REAL_FUZZ);

    VecF samples(10, 0);

    const int N = 10000000;
    for (int i = 0; i < N; i++) {
        samples[rng.sample_cdf(cdf)]++;
    }

    std::ofstream sample_out("pdf.txt");
    for (unsigned i = 0; i < pdf.size(); i++) {
        std::cout << samples[i] / N << std::endl;
        sample_out << samples[i] / N << " " << pdf[i] << std::endl;
    }
}

int main()
{
    return UnitTest::RunAllTests();
}
