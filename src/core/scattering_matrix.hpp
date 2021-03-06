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

#include <iosfwd>
#include "util/blitz_typedefs.hpp"
#include "util/global_config.hpp"

namespace mocc {
// Scattering matrix row
struct ScatteringRow {
public:
    ScatteringRow(int min, int max, real_t const *const from)
        : min_g(min), max_g(max), from(from)
    {
        return;
    }

    int min_g;
    int max_g;
    real_t const *const from;

    real_t operator[](int g) const
    {
        assert(g >= min_g);
        assert(g <= max_g);
        return from[g - min_g];
    }

    const real_t *begin() const
    {
        return from;
    }

    const real_t *end() const
    {
        return from + max_g - min_g + 1;
    }

    bool operator==(const ScatteringRow &other) const
    {
        if (this == &other) {
            return true;
        }
        if (min_g != other.min_g) {
            return false;
        }
        if (max_g != other.max_g) {
            return false;
        }

        return true;
    }

    bool operator!=(const ScatteringRow &other) const
    {
        return !(*this == other);
    }
};

/**
 * This class provides an efficient means by which to store a matrix of
 * scattering cross sections. Generally speaking, scattering matrices tend
 * to be relatively sparse, since upscatter is not present at high energies
 * (so the matrix is largely lower-triangular), and downscattering energy
 * transfer is physically limited by the ratio of masses. Therefore we use a
 * compressed representation, where each "row" of outscatter cross sections
 * are stored contiguously, along with their group boundaries.
 */
class ScatteringMatrix {
public:
    ScatteringMatrix() : ng_(0)
    {
    }

    /**
     * Construct a scattering matrix using a 2-dimensional vector<vector<>>.
     * This full, dense representation of the scattering matrix will be
     * sparsified.
     *
     * \param scat the dense representation of the scattering matrix.
     * Indexing should be [to group][from group]
     */
    ScatteringMatrix(const std::vector<VecF> &scat);

    /**
     * Construct a scattering matrix using a 2-dimensional Blitz++ array.
     * This full, dense representation of the scattering matrix will be
     * sparsified.
     *
     * \param scat the dense representation of the scattering matrix.
     * Indexing should be (to group, from group)
     */
    ScatteringMatrix(const ArrayB2 &scat);

    /**
     * Copy constructor. Need this in order to produce valid raw pointers to
     * the scattering rows.
     */
    ScatteringMatrix(const ScatteringMatrix &other)
        : ng_(other.ng_), scat_(other.scat_), out_(other.out_)
    {
        // Pretty much everything can copy straight over, but we need to
        // reach into the scattering rows and update their pointers to the
        // location of the new scat_ vector
        int pos = 0;
        for (auto &row : other) {
            rows_.push_back(ScatteringRow(row.min_g, row.max_g, &scat_[pos]));
            pos += row.max_g - row.min_g + 1;
        }

        return;
    }

    ScatteringMatrix &operator=(const ScatteringMatrix &rhs)
    {
        if (this != &rhs) {
            ng_   = rhs.ng_;
            scat_ = rhs.scat_;
            out_  = rhs.out_;
            rows_.clear();

            int pos = 0;
            for (auto &row : rhs) {
                rows_.push_back(
                    ScatteringRow(row.min_g, row.max_g, &scat_[pos]));
                pos += row.max_g - row.min_g + 1;
            }
        }
        return *this;
    }

    const ScatteringRow &to(int ig) const
    {
        assert((ig >= 0) && (ig < int(rows_.size())));
        return rows_[ig];
    }

    /**
     * \brief Return the self-scattering cross section for the indicated
     * group.
     */
    real_t self_scat(int group) const
    {
        return rows_[group][group];
    }

    /**
     * Return the number of energy groups for which the scattering matrix is
     * defined.
     */
    int n_group() const
    {
        return ng_;
    }

    /**
     * Return the total out-scattering cross section for group ig
     *
     * This includes self-scatter, and is equivalent to a column sum of the
     * full scattering matrix.
     *
     * \todo this is improperly-named. Should just be "total" or something
     */
    real_t out(unsigned int ig) const
    {
        return out_[ig];
    };

    /**
     * \brief Return a CDF of the outscatter probabilities for group \p ig
     */
    VecF out_cdf(int ig) const
    {
        VecF cdf;
        cdf.reserve(ng_);

        real_t scale = 1.0 / this->out(ig);
        real_t prev  = 0.0;
        for (int igg = 0; igg < ng_; igg++) {
            real_t xssc =
                ((ig >= this->to(igg).min_g) && (ig <= this->to(igg).max_g))
                    ? this->to(igg)[ig]
                    : 0.0;
            cdf.push_back(prev + xssc * scale);
            prev = cdf.back();
        }
        return cdf;
    }

    /**
     * Return iterator to the first scattering row.
     */
    std::vector<ScatteringRow>::const_iterator begin() const
    {
        return rows_.cbegin();
    }

    /**
     * Return iterator past the last scattering row.
     */
    std::vector<ScatteringRow>::const_iterator end() const
    {
        return rows_.cend();
    }

    /**
     * \brief Return a 1-D, dense representation of the scattering matrix.
     *
     * The returned vector stores all of scattering cross sections as a
     * row-major ng-by-ng matrix.
     */
    VecF as_vector() const
    {
        VecF sc(ng_ * ng_, 0.0);
        int ig = 0;
        for (const auto &row : rows_) {
            for (int igg = row.min_g; igg <= row.max_g; igg++) {
                sc[ng_ * ig + igg] = row[igg];
            }
            ig++;
        }
        return sc;
    }

    bool operator==(const ScatteringMatrix &other) const
    {
        if (this == &other) {
            return true;
        }
        if (ng_ != other.n_group()) {
            return false;
        }
        if (scat_ != other.scat_) {
            return false;
        }
        if (out_ != other.out_) {
            return false;
        }
        return true;
    }

    bool operator!=(const ScatteringMatrix &other) const
    {
        return !(*this == other);
    }

    // Provide stream insertion support
    friend std::ostream &operator<<(std::ostream &os,
                                    const ScatteringMatrix &scat_mat);

private:
    int ng_;
    // Densified scattering cross sections
    VecF scat_;
    // Group-wise outscatter cross sections
    VecF out_;
    std::vector<ScatteringRow> rows_;
};
}
