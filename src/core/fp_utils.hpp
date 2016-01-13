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

    inline bool fp_equiv_ulp(float v1, float v2) {
        fp_utils::float_int i1;
        i1.f = v1;
        if (i1.i < 0) {
            i1.i = 0x80000000 - i1.i;
        }

        fp_utils::float_int i2;
        i2.f = v2;
        if (i2.i < 0) {
            i2.i = 0x80000000 - i2.i;
        }

        return std::abs(i1.i - i2.i) < 100;
    }

    inline bool fp_equiv_ulp(double v1, double v2) {
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

        return std::abs(i1.i - i2.i) < 100;
    }

    inline bool fp_equiv_rel(real_t v1, real_t v2) {
        return fabs(v1-v2)/fabs(v1) < FLOAT_EPS;
    }
    inline bool fp_equiv_abs(real_t v1, real_t v2) {
        return fabs(v1-v2) < FLOAT_EPS;
    }
}
