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

#include <cassert>

#include "global_config.hpp"
#include "util/force_inline.hpp"

namespace mocc {
    /**
     * \brief A linear congruential random number generator
     *
     * Most of the parameters are from OpenMC. Thanks, dudes!
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

        /**
         * \brief Generate a uniformly-distributed random number on
         * [0,\p ubound)
         */
        MOCC_FORCE_INLINE real_t random( real_t ubound ) {
            assert( ubound > 0.0 );
            current_seed_ = current_seed_*m_ + increment_;

            real_t v = float_scale_ * current_seed_;
            return v*ubound;
        }

        /**
         * \brief Generate a uniformly-distributed random number on
         * [\p lbound, \p ubound)
         */
        MOCC_FORCE_INLINE real_t random( real_t lbound, real_t ubound ) {
            assert( ubound > lbound );
            current_seed_ = current_seed_*m_ + increment_;

            real_t v = float_scale_ * current_seed_;
            return lbound + (ubound-lbound)*v;
        }

        /**
         * \brief Sample a uniformly-distributed integer on [0, ubound)
         */
        MOCC_FORCE_INLINE int random_int( int ubound ) {
            current_seed_ = current_seed_*m_ + increment_;
            int i = current_seed_ * ubound;

            return i;
        }

        /**
         * \brief Sample and index from a cumulative distribution function
         *
         * \param cdf a vector of \ref real_t containing the CDF. This is
         * assumed to increase monotonically to a final value of unity.
         *
         * This will sample an index randomly from a CDF, which should contain
         * the cumulative probability of the index lying below each entry. To be
         * well-formed, the entries in the CDF should increase monotonically,
         * with the last entry being unity. There is no check made internally,
         * and it is assuming that the caller is providing a valid CDF.
         */
        MOCC_FORCE_INLINE int sample_cdf( const std::vector<real_t> &cdf ) {
            real_t v = this->random();
            return std::distance(cdf.begin(),
                    std::lower_bound(cdf.begin(), cdf.end(), v) );
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

    extern RNG_LCG RNG;

} // namespace mocc
