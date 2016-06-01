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
 * \brief Monte Carlo tally for a spatially-dependent quantity
 *
 * See \ref tally_page for more discussion about tallies
 */
class TallySpatial {
public:
    /**
     * \brief Make a new \ref TallyScalar
     */
    TallySpatial(const VecF &norm)
        : nreg_(norm.size()), norm_(norm), data_(nreg_, {0.0, 0.0}), weight_(0.0)
    {
        return;
    }

    /**
     * \brief Score some quantity to the tally
     */
    void score(int i, real_t value)
    {
#pragma omp atomic
        data_[i].first += value;
#pragma omp atomic
        data_[i].second += value*value;

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
        for (auto &d : data_) {
            d = {0.0, 0.0};
        }
        weight_ = 0.0;
        return;
    }

    /**
     * \brief Return the estimates for the tally mean and variance
     */
    const std::vector<std::pair<real_t, real_t>> get() const
    {
        std::vector<std::pair<real_t, real_t>> ret;
        ret.reserve(data_.size());
        for( unsigned i=0; i<data_.size(); i++) {
            real_t mean           = data_[i].first / (weight_*norm_[i]);
            real_t mean_of_square = data_[i].second / (weight_*norm_[i]);
            real_t square_of_mean = mean * mean;
            real_t variance =
                (mean_of_square - square_of_mean) / (weight_ - 1.0);

            ret.push_back({mean, variance});
        }

        return ret;
    }

private:
    unsigned nreg_;
    const VecF &norm_;
    std::vector<std::pair<real_t, real_t>> data_;
    real_t weight_;
};
}
} // namespaces
