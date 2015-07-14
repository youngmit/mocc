#pragma once

#include "global_config.hpp"

namespace mocc {

    inline bool fp_equiv_ulp(float_t v1, float_t v2) {
        int i1 = *(int*) &v1;
        if (i1 < 0) {
            i1 = 0x80000000 - i1;
        }
        int i2 = *(int*) &v2;
        if (i2 < 0) {
            i2 = 0x80000000 - i2;
        }

        return abs(i1 - i2) < 100; 
    }

    inline bool fp_equiv_rel(float_t v1, float_t v2) {
        return fabs(v1-v2)/fabs(v1) < FLOAT_EPS;
    }
    inline bool fp_equiv_abs(float_t v1, float_t v2) {
        return fabs(v1-v2) < FLOAT_EPS;
    }
}
