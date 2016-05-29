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

#include "util/rng_lcg.hpp"

namespace mocc {
namespace mc {
class RNGSwarm {
public:
    RNGSwarm(int rank, int width, unsigned long int seed = 1)
        : master_(seed), generators_(), rank_(rank), width_(width)
    {
        assert(rank >= 0);
        assert(width > 0);
        // Advance the master generator up by rank*width;
        master_.jump_ahead(rand_ * width_);
        generators_.reserve(width_);
        for (int i = 0; i < width_; i++) {
            generators_.emplace_back(master_.get());
        }
        return;
    }

    RNG_LCG &operator[](int i){
        return generators_[i];
    }

private:
    RNG_LCG master_;
    std::vector<RNG_LCG> generators_;
};
}
}
