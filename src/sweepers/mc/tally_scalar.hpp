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
    void score(real_t tl, real_t multiplier)
    {
        real_t contrib = tl * multiplier;
        mean_ += contrib;
        mean_square_ += contrib * contrib;

        return;
    }

    /**
     * \brief Introduce new weight to the tally
     */
    void add_weight(real_t w)
    {
        weight_ += w;
    }

    /**
     * \breif Reset the tally, forgetting all history
     */
    void reset()
    {
        mean_        = 0.0;
        mean_square_ = 0.0;
        weight_      = 0.0;
        return;
    }

    /**
     * \brief Return the estimates for the tally mean and variance
     */
    std::pair<real_t, real_t> get() const
    {
std::cout << mean_ << " " << mean_square_ << " " << weight_ << std::endl;
        std::pair<real_t, real_t> val;
        real_t mean = mean_ / weight_;
        real_t mean_of_square = mean_square_/(weight_);
        real_t square_of_mean = mean * mean;
        real_t variance = (mean_of_square - square_of_mean)/(weight_ - 1.0 );

        val.first = mean;
        val.second = variance;
        return val;
    }

private:
    real_t mean_;
    real_t mean_square_;
    real_t weight_;
};
}
} // namespaces
