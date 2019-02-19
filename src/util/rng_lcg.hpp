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
#include <cstdint>
#include <iostream>
#include <random>

#include "util/force_inline.hpp"
#include "global_config.hpp"
#include "fp_utils.hpp"

namespace mocc {
/**
 * \brief A linear congruential random number generator
 *
 * Most of the parameters are from OpenMC/MCNP random number generators
 */
class RNG_LCG {
public:
    RNG_LCG(uint64_t seed = UINT64_C(1)) : current_seed_(seed)
    {
        return;
    }

    RNG_LCG &operator=(const RNG_LCG &other)
    {
        // Skip self-assignment check, since we only have POD data members
        current_seed_ = other.current_seed_;
        return *this;
    }

    void set_seed(uint64_t seed)
    {
        current_seed_ = seed;
    }

    uint64_t operator()()
    {
        current_seed_ = (current_seed_ * m_ + b_) & mask_;
        return current_seed_;
    }

    /**
     * \brief Generate a uniformly-distributed random number on [0,1)
     */
    real_t random()
    {
        return float_scale_ * (*this)();
    }

    /**
     * \brief Generate a uniformly-distributed random number on
     * [0,\p ubound)
     */
    MOCC_FORCE_INLINE real_t random(real_t ubound)
    {
        assert(ubound > 0.0);

        real_t v = float_scale_ * (*this)();
        return v * ubound;
    }

    /**
     * \brief Generate a uniformly-distributed random number on
     * [\p lbound, \p ubound)
     */
    MOCC_FORCE_INLINE real_t random(real_t lbound, real_t ubound)
    {
        assert(ubound > lbound);

        real_t v = this->random();
        v        = lbound + (ubound - lbound) * v;
        return v;
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
        assert(fp_equiv_ulp(cdf.back(), 1.0));
        real_t v = this->random();
        return std::distance(cdf.begin(),
                             std::lower_bound(cdf.begin(), cdf.end(), v));
    }

    /**
     * \brief Move the state of the generator forward in the sequence
     *
     * \param n the number of elements in the sequence to jump ahead by
     *
     */
    void jump_ahead(uint64_t n)
    {
        uint64_t nskip = n;

        while (nskip < 0l) {
            nskip += mod_;
        }

        uint64_t g     = m_;
        uint64_t b     = b_;
        uint64_t g_new = 1;
        uint64_t b_new = 0;
        while (nskip > UINT64_C(0)) {
            if (nskip & UINT64_C(1)) {
                g_new = g_new * g & mask_;
                b_new = (b_new * g + b) & mask_;
            }

            b = ((g + 1) * b) & mask_;
            g = (g * g) & mask_;

            // Shift bits left
            nskip = nskip >> 1;
        }

        current_seed_ = (g_new * current_seed_ + b_new) & mask_;
        return;
    }

private:
    uint64_t current_seed_;
    static const uint64_t m_ = UINT64_C(2806196910506780709);
    static const uint64_t b_ = 1;
    // The mask performs the equivalent of a modulo (%) when ANDed with
    // left-hand operand.
    static const uint64_t mask_     = ~(UINT64_C(1) << 63);
    static const uint64_t mod_      = UINT64_C(1) << 63;
    static constexpr real_t float_scale_ = 1.0 / (UINT64_C(1) << 63);
};

} // namespace mocc
