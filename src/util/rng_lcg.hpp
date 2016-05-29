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
#include <iostream>
#include <random>

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
    RNG_LCG(unsigned long seed = 1ul) : generator_(seed)
    {
        return;
    }

    RNG_LCG &operator=(const RNG_LCG &other)
    {
        // Skip self-assignment check, since we only have POD data members
        generator_ = other.generator_;
        return *this;
    }

    unsigned long operator()()
    {
        return generator_();
    }

    /**
     * \brief Generate a uniformly-distributed random number on [0,1)
     */
    real_t random()
    {
        return float_scale_ * generator_();
    }

    /**
     * \brief Generate a uniformly-distributed random number on
     * [0,\p ubound)
     */
    MOCC_FORCE_INLINE real_t random(real_t ubound)
    {
        assert(ubound > 0.0);

        real_t v = float_scale_ * generator_();
        return v * ubound;
    }

    /**
     * \brief Generate a uniformly-distributed random number on
     * [\p lbound, \p ubound)
     */
    MOCC_FORCE_INLINE real_t random(real_t lbound, real_t ubound)
    {
        assert(ubound > lbound);

        real_t v = float_scale_ * generator_();
        assert(v >= 0.0);
        assert(v < 1.0);
        return lbound + (ubound - lbound) * v;
    }

    /**
     * \brief Sample a uniformly-distributed integer on [0, ubound)
     */
    MOCC_FORCE_INLINE int random_int(int ubound)
    {
        int i = this->random() * ubound;

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
    MOCC_FORCE_INLINE int sample_cdf(const std::vector<real_t> &cdf)
    {
        real_t v = this->random();
        return std::distance(cdf.begin(),
                             std::lower_bound(cdf.begin(), cdf.end(), v));
    }

    /**
     * \brief Move the state of the generator forward in the sequence
     *
     * \param n the number of elements in the sequence to jump ahead by
     *
     * This is using the stupid approach for now. \todo use the logarithmic
     * approach
     */
    void jump_ahead(int n)
    {
        for (int i = 0; i < n; i++) {
            generator_();
        }
    }

private:
    std::linear_congruential_engine<unsigned long, 2806196910506780709ul, 1ul,
                                    0>
        generator_;
    static constexpr real_t float_scale_ = 1.0 /
        std::linear_congruential_engine<unsigned long, 2806196910506780709ul,
                                        1ul, 0>::max();
};

} // namespace mocc
