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

#include "core/global_config.hpp"

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
    TallyScalar() : mean_(0.0), mean_square_(0.0), weight_(0.0)
    {
        return;
    }

    /**
     * \brief Score some quantity to the tally
     */
    void score(real_t value)
    {
#pragma omp atomic
        mean_ += value;
#pragma omp atomic
        mean_square_ += value * value;

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
            mean_        = 0.0;
            mean_square_ = 0.0;
            weight_      = 0.0;
        }

        return;
    }

    /**
     * \brief Return the estimates for the tally mean and relative standard
     * deviation
     */
    std::pair<real_t, real_t> get() const
    {
        // Plenty of room for optimization in here
        std::pair<real_t, real_t> val;
        real_t mean           = mean_ / weight_;
        real_t mean_of_square = mean_square_ / (weight_ - 1.0);
        real_t square_of_mean = mean_ * mean_ / weight_ * (weight_ - 1.0);
        real_t variance = mean_of_square - square_of_mean;

        val.first  = mean;
        val.second = std::sqrt(variance) / mean;
        return val;
    }

private:
    real_t mean_;
    real_t mean_square_;
    real_t weight_;
};
}
} // namespaces
