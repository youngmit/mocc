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

using std::cout;
using std::endl;

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
     * \brief Score some quantity to the tally
     *
     * This contributes to the running sum for the mean and the sum of squares.
     * This is correct for collision-style estimates, but yeilds incorrect
     * statistical estimates for track length-based elements, since the tallies
     * for each entire history should be squared, not each individual track. For
     * now we use batch statistics to get around this, but at some point it
     * might be nice to do track length tallies directly.
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
            real_t r_weight = 1.0/weight_;
            for (unsigned i = 0; i < realization_scores_.size(); i++) {
                real_t v = realization_scores_[i] * r_weight;
                data_[i].first += v;
                data_[i].second += v*v;
                realization_scores_[i] = 0.0;
            }
            n_++;
            weight_ = 0.0;
        }
    }

    /**
     * \brief Score the contents of an entire tally to this tally.
     *
     * This is particularly useful for computing batch statistics
     */
    void score(const TallySpatial &other)
    {
#pragma omp single
        {

throw EXCEPT("shouldnt be using this anymore");
            assert(data_.size() == other.data_.size());
            real_t w = other.weight_;
            for (unsigned i = 0; i < data_.size(); i++) {
                real_t v = other.data_[i].first / (w);
                data_[i].first += v;
                data_[i].second += v * v;
            }
        }
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
            for (auto &d : data_) {
                d = {0.0, 0.0};
            }
            for (auto &d : realization_scores_) {
                d = 0.0;
            }
            weight_ = 0.0;
            n_ = 0;
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
        assert(weight_ > 0.0);
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
