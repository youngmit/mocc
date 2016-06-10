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

#pragma once

#include <math.h>
#include <vector>

#include "util/error.hpp"
#include "util/global_config.hpp"
#include "core/geometry/angle.hpp"

using namespace mocc;

// Hide all of these constants from the rest of the world
namespace {
real_t mu_base[] = {0.577350269189626, 0.350021000000000, 0.266636000000000,
                    0.218218218218218, 0.192450089729876, 0.174077655955702,
                    0.161575000000000, 0.149071198499989};

real_t w_unique[] = {
    1.0,       1.0 / 3.0, 0.1761262, 0.1572071, 0.1209876, 0.0907408, 0.0925925,
    0.0893043, 0.0725281, 0.0450455, 0.0539274, 0.0707734, 0.0558760, 0.0373436,
    0.0502654, 0.0258553, 0.0580031, 0.0488943, 0.0228095, 0.0393955, 0.0380920,
    0.0258382, 0.0082759, 0.0489967, 0.0413235, 0.0203158, 0.0265468, 0.0378883,
    0.0135404, 0.0326129, 0.0103825};

int w_offset[] = {0, 1, 2, 4, 7, 11, 16, 23};

int w_map[] = {1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 2, 2, 3, 2, 1, 2, 2, 1,
               1, 2, 2, 3, 4, 3, 2, 4, 4, 2, 1, 2, 3, 2, 1, 1, 2, 2, 3, 5,
               3, 4, 6, 6, 4, 3, 6, 7, 6, 3, 2, 5, 6, 6, 5, 2, 1, 2, 3, 4,
               3, 2, 1, 1, 2, 2, 3, 5, 3, 4, 6, 6, 4, 4, 7, 8, 7, 4, 3, 6,
               8, 8, 6, 3, 2, 5, 6, 7, 6, 5, 2, 1, 2, 3, 4, 4, 3, 2, 1};

int w_map_offset[] = {0, 1, 4, 10, 20, 35, 63};
}

// Produce a vector of angles matching the level-symmetric quadrature of order
// 'order'
std::vector<Angle> GenSn(int order)
{
    if (order % 2 != 0) {
        throw EXCEPT("Sn quadrature order must be even.");
    }
    if (order > 16) {
        throw EXCEPT("Max supported Sn quadrature order is 16.");
    }

    // n is the number of base cosines to use
    const int n = order / 2;

    // set up the list of base cosines
    VecF mu;
    mu.push_back(mu_base[n - 1]);
    if (order > 2) {
        const real_t delta_mu =
            2.0 * (1.0 - 3.0 * (mu[0] * mu[0])) / (real_t)(order - 2);
        for (int i = 1; i < n; i++) {
            mu.push_back(sqrt(mu[0] * mu[0] + i * delta_mu));
        }
    }

    // Alias the w_unique array to get a slice for the order we are interested
    // in.
    real_t *weights = &w_unique[w_offset[n - 1]];
    // Alias into the w_map to get our indices
    int *map = &w_map[w_map_offset[n - 1]];

    // Apply the permutations of the base cosines to get actual angles. We will
    // do this once for the first octant, then reflect around.
    std::vector<Angle> angles;
    real_t wsum = 0.0;
    int k       = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            Angle angle(mu[i - j], mu[j], mu[n - i - 1], weights[map[k] - 1]);

            angles.push_back(angle);
            wsum += angle.weight;
            k++;
        }
    }

    // One of these days, i should calculate the weights algorithmically, rather
    // than storing in a table. For now, make sure the angular integral comes
    // out to 4*PI. We defined the angles above for one octant, so make sure
    // that the weights sum to 1 to machine precision.
    for (auto &a : angles) {
        a.weight /= wsum;
    }

    return angles;
}
