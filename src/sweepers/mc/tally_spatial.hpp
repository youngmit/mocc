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
 * \brief Monte Carlo tally for a spatially-dependent quantity
 *
 * Calls to \ref score() contribute to a buffer, \ref realization_scores_, which
 * following the completion of a "sample" can then stored to the persistent
 * tally values, \ref data_, using the \ref commit_realization() method. \ref
 * data_ stores a sequence of \c std::pair, each containing a running sum and
 * sum of the square of the values from each realization for a region of phase
 * space.
 *
 * Calling \ref get() returns the mean and relative standard deviation for each
 * region of phase space.
 *
 * See \ref tally_page for more discussion about tallies
 */
class TallySpatial {
public:
    /**
     * \brief Make a new \ref TallyScalar
     */
    TallySpatial(const VecF &norm)
        : nreg_(norm.size()),
          norm_(norm),
          data_(nreg_, {0.0, 0.0}),
          realization_scores_(norm.size(), 0.0),
          weight_(0.0),
          n_(0)
    {
        return;
    }

    /**
     * \brief Score some quantity to the tally realization buffer
     */
    void score(int i, real_t value)
    {
#pragma omp atomic
        realization_scores_[i] += value;

        return;
    }

    /**
     * \brief Commit tally contributions for a given realization to the tally
     */
    void commit_realization()
    {
#pragma omp single
        {
            real_t r_weight = 1.0 / weight_;
            for (unsigned i = 0; i < realization_scores_.size(); i++) {
                real_t v = realization_scores_[i] * r_weight;
                data_[i].first += v;
                data_[i].second += v * v;
                realization_scores_[i] = 0.0;
            }
            n_++;
            weight_ = 0.0;
        }
    }

    /**
     * \brief Introduce new weight to the tally
     */
    void add_weight(real_t w)
    {
#pragma omp atomic
        weight_ += w;
        return;
    }

    /**
     * \brief Reset the tally, forgetting all history
     */
    void reset()
    {
#pragma omp single
        {
            for (auto &d : data_) {
                d = {0.0, 0.0};
            }
            for (auto &d : realization_scores_) {
                d = 0.0;
            }
            weight_ = 0.0;
            n_      = 0;
        }
        return;
    }

    /**
     * \brief Return the estimates for the tally mean and relative standard
     * deviation
     */
    const std::vector<std::pair<real_t, real_t>> get() const
    {
        // Plenty of room for optimization in here
        assert(n_ > 0);
        std::vector<std::pair<real_t, real_t>> ret;
        ret.reserve(data_.size());
        for (unsigned i = 0; i < data_.size(); i++) {
            real_t r_norm = 1.0 / norm_[i];
            real_t mean   = data_[i].first / n_;
            real_t square_of_mean =
                data_[i].first * data_[i].first / (n_ * (n_ - 1.0));
            real_t mean_of_square = data_[i].second / (n_ - 1.0);
            real_t variance       = mean_of_square - square_of_mean;

            ret.push_back({mean * r_norm,
                           std::sqrt(variance * r_norm) / (mean * r_norm)});
        }

        return ret;
    }

private:
    unsigned nreg_;
    const VecF &norm_;
    std::vector<std::pair<real_t, real_t>> data_;
    VecF realization_scores_;
    real_t weight_;
    int n_;
};
}
} // namespaces
