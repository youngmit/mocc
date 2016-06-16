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

#include <iostream>
#include <utility>
#include "util/global_config.hpp"

namespace mocc {
namespace mc {

/**
 * \brief Monte Carlo tally for a scalar quantity
 *
 * See \ref tally_page for more discussion about tallies
 */
class TallyScalar {
public:
    /**
     * \brief Make a new \ref TallyScalar
     */
    TallyScalar() : sum_(0.0), sum_square_(0.0), weight_(0.0)
    {
        return;
    }

    /**
     * \brief Score some quantity to the tally
     */
    void score(real_t value)
    {
#pragma omp atomic
        sum_ += value;
#pragma omp atomic
        sum_square_ += value * value;

        return;
    }

    /**
     * \brief Introduce new weight to the tally
     */
    void add_weight(real_t w)
    {
#pragma omp atomic
        weight_ += w;
    }

    /**
     * \brief Reset the tally, forgetting all history
     */
    void reset()
    {
#pragma omp single
        {
            sum_        = 0.0;
            sum_square_ = 0.0;
            weight_     = 0.0;
        }

        return;
    }

    /**
     * \brief Return the estimates for the tally mean and relative standard
     * deviation of the mean
     */
    std::pair<real_t, real_t> get() const
    {
        // Plenty of room for optimization in here
        std::pair<real_t, real_t> val;
        real_t mean           = sum_ / weight_;

        val.first  = mean;
        val.second = 1.0 / (weight_ - 1.0) *
                     ((1.0 / weight_) * sum_square_ - mean * mean);
        val.second = std::sqrt(val.second) / mean;
        return val;
    }

private:
    real_t sum_;
    real_t sum_square_;
    real_t weight_;
};
}
} // namespaces
