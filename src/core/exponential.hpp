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

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "util/global_config.hpp"

namespace mocc {
/**
 * Base class defining a class of simple utility classes for computing
 * exponential functions. The base version uses the stock library call to
 * exp(), while derived versions can override this with more efficient table
 * lookups.
 */
class Exponential {
public:
    Exponential()
    {
    }

    inline real_t exp(real_t v) const
    {
        return std::exp(v);
    }

    virtual real_t max_error()
    {
        return 0.0;
    }
};

/*
 * This version of \ref Exponential uses a linear lookup table to speed up
 * the evaluations of \ref exp(). If the argument to \ref exp() is beyond
 * the domain of the table, it will fall back to the result of the standard
 * library \c exp() function.
 */
template <int N> class Exponential_Linear : public Exponential {
public:
    Exponential_Linear(real_t min = -10.0, real_t max = 0.0)
        : min_(min),
          max_(max),
          space_((max - min_) / (real_t)(N)),
          rspace_(1.0 / space_)
    {
        for (int i = 0; i <= N; i++) {
            d_[i] = std::exp(min_ + i * space_);
        }
    }

    inline real_t exp(real_t v) const
    {
        if (v < min_ || v > max_) {
            std::cout << "Out-of-bounds exponential argument: " << v
                      << std::endl;
            return std::exp(v);
        }
        int i = (v - min_) * rspace_;
        v -= space_ * i + min_;
        return d_[i] + (d_[i + 1] - d_[i]) * v * rspace_;
    }

    real_t max_error()
    {
        real_t max_error = 0.0;
        for (int i = 0; i < N; i++) {
            real_t x   = min_ + space_ * (0.5 + i);
            real_t e   = std::exp(x);
            real_t err = std::abs((this->exp(x) - e) / e);
            max_error  = std::max(max_error, err);
        }
        return max_error;
    }

    /*
     * \brief Return the table value for the passed point index
     *
     * This is mostly useful for testing and debugging purposes.
     */
    real_t operator[](int i) const {
        return d_[i];
    }

    real_t dx() const {
        return space_;
    }

protected:
    real_t min_;
    real_t max_;
    real_t space_;
    real_t rspace_;
    std::array<real_t, N + 1> d_;
};

/**
 * Same as \ref Exponential_Linear, but without the check to make sure the
 * argument is within the limits of the table. This should only be used in
 * situations where one knows that the arguments to \ref exp() won't spill
 * the banks of the table. Even considering branch prediction, this manages
 * to shave a little more time off of the \ref exp() evaluations over \ref
 * Exponential_Linear.
 */
template <int N> class Exponential_UnsafeLinear : public Exponential_Linear<N> {
    Exponential_UnsafeLinear(real_t min = -10.0, real_t max = 0.0)
        : Exponential_Linear<N>(min, max)
    {
        return;
    }

    inline real_t exp(real_t v) const
    {
        int i = (v - this->min_) * this->rspace_;
        v -= this->space_ * i + this->min_;
        return this->d_[i] +
               (this->d_[i + 1] - this->d_[i]) * v * this->rspace_;
    }
};
}
