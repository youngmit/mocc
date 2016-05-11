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

#include "global_config.hpp"
#include "util/force_inline.hpp"

namespace mocc {
    /**
     * \brief A linear congruential random number generator
     */
    class RNG_LCG {
    public:
        RNG_LCG( unsigned long seed=1 ):
            seed_(seed)
        {
            return;
        }

        /**
         * \brief Generate a uniformly-distributed random number on [0,1)
         */
        MOCC_FORCE_INLINE real_t random() {
            current_seed_ = current_seed_*m_ + increment_;

            return float_scale_ * current_seed_;
        }
        
    private:
        const unsigned long seed_;
        unsigned long current_seed_;
        const int bits_ = 64;
        const unsigned long m_ = 2806196910506780709ul; 
        const unsigned long mod_ = -1;
        const unsigned long increment_ = 1;
        const real_t float_scale_ = 1.0/(std::pow(2.0, bits_));
    };
} // namespace mocc
