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
#include <cmath>
#include "global_config.hpp"

namespace mocc {
namespace fp_utils {
union float_int {
    float f;
    int32_t i;
};

union double_int {
    double f;
    int64_t i;
};
}

/// \todo this kind of thing is multiply-defined throughout the code. Try
/// and settle on one value
constexpr real_t REAL_FUZZ = 10.0 * std::numeric_limits<real_t>::epsilon();

/**
 * \brief Compare two floats using ULP.
 *
 * See <a
 * href="http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm">
 * this page</a> and <a
 * href="https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
 * this page</a> for more info about testing floating point equivalence.
 *
 * \note This routine is safe for comparing values that are not terribly
 * close to zero. For comparing near zero, use \ref fp_equiv_abs().
 */
inline bool fp_equiv_ulp(real_t v1, real_t v2)
{
    fp_utils::double_int i1;
    i1.f = v1;
    if (i1.i < 0) {
        i1.i = 0x80000000 - i1.i;
    }

    fp_utils::double_int i2;
    i2.f = v2;
    if (i2.i < 0) {
        i2.i = 0x80000000 - i2.i;
    }

    return std::abs(i1.i - i2.i) < 200;
}

inline bool fp_equiv_rel(real_t v1, real_t v2)
{
    return std::abs(v1 - v2) / std::abs(v1) < REAL_FUZZ;
}
inline bool fp_equiv_abs(real_t v1, real_t v2)
{
    return std::abs(v1 - v2) < REAL_FUZZ;
}

/**
 * \brief Compare two floating point values for aproximate equivalence.
 *
 * This one is the kitchen sink; if the two numbers are close by absolute
 * comparison, return true, otherwise apply a ULP-based comparison. This is one
 * of the methods suggested by randomascii.
 */
inline bool fp_equiv(real_t v1, real_t v2)
{
    if (fp_equiv_abs(v1, v2)) {
        return true;
    }

    // If the signs differ, not equal
    if ((v1 < 0.0) != (v2 < 0.0)) {
        return false;
    }

    return fp_equiv_ulp(v1, v2);
}

/**
 * \brief A lambda function for doing a fuzzy floating-point less-than (\<)
 * comparison.
 *
 * This function returns true when \p l is sufficiently less than \p r. The
 * primary utility of such a function is for use with the \c
 * std::lower_bound() and \c std::upper_bound() algorithms.
 */
auto fuzzy_lt = [](real_t l, real_t r) { return (l - r) < -REAL_FUZZ; };
}
