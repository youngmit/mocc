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

    /**
     * \brief Compare two floats using ULP.
     *
     * See <a
     * href="http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm">
     * this page</a> and <a
     * href="https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/">
     * this page</a> for more info about testing floating point equivalence.
     */
    inline bool fp_equiv_ulp(real_t v1, real_t v2) {
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

    inline bool fp_equiv_rel(real_t v1, real_t v2) {
        return fabs(v1-v2)/fabs(v1) < FLOAT_EPS;
    }
    inline bool fp_equiv_abs(real_t v1, real_t v2) {
        return fabs(v1-v2) < FLOAT_EPS;
    }
}
