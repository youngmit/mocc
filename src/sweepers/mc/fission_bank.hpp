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

#include <array>
#include <iosfwd>
#include <vector>

#include "util/global_config.hpp"
#include "util/pugifwd.hpp"
#include "util/rng_lcg.hpp"
#include "core/core_mesh.hpp"
#include "core/geometry/geom.hpp"
#include "core/xs_mesh.hpp"
#include "particle.hpp"

namespace mocc {
namespace mc {
/**
 * A FissionBank stores a sequence of fission sites. Nothing fancy
 */
class FissionBank {
public:
    FissionBank(const CoreMesh &mesh);

    /**
     * \brief Construct a FissionBank by uniformly sampling fission sites.
     *
     * \param input XML node containing bounds of a 3-D box within which to
     * sample initial fission sites
     * \param n the number of initial sites to sample
     * \param mesh the \ref CoreMesh to use for initial sampling
     * \param xs_mesh the \ref XSMesh to use for initial sampling
     * \param rng a reference to the random number generator to be used for
     * sampling initial fission sites.
     *
     * This constructor initializes a \ref FissionBank using input specified in
     * an XML node.
     */
    FissionBank(const pugi::xml_node &input, int n, const CoreMesh &mesh,
                const XSMesh &xs_mesh, RNG_LCG &rng);

    auto begin()
    {
        return sites_.begin();
    }

    const auto begin() const
    {
        return sites_.cbegin();
    }

    auto end()
    {
        return sites_.end();
    }

    const auto end() const
    {
        return sites_.cend();
    }

    int size() const
    {
        return sites_.size();
    }

    /**
     * \brief Add a new fission site to the \ref FissionBank
     *
     * \param p a \ref Point3 for the location of the fission site
     *
     * This method adds a new fission site to the fission bank, and makes a
     * contribution to the total number of neutrons that were generated into
     * the bank.
     */
    void push_back(Particle &p)
    {
#pragma omp critical
        {
            sites_.push_back(p);
            total_fission_ += p.weight;
        }
        return;
    }

    /**
     * \brief Return the Shannon entropy of the fission bank.
     *
     * This is used to estimate the change in the spatial distribution of
     * fission sites from generation to generation. Observing little
     * variation in this metric throughout the active cycles lends some
     * confidence that the fission source distribution was well converged
     * before beginning active cycles.
     */
    real_t shannon_entropy() const;

    /**
     * \brief Swap contents with another \ref FissionBank
     *
     * \param other the other \ref FissionBank to swap with
     */
    void swap(FissionBank &other);

    const auto &operator[](unsigned i) const
    {
        return sites_[i];
    }

    /**
     * \brief Clear the \ref FissionBank of all fission sites
     */
    void clear()
    {
#pragma omp single
        {
            sites_.clear();
            total_fission_ = 0.0;
        }
    }

    void resize(unsigned int n, RNG_LCG &rng);

    real_t total_fission() const
    {
        return total_fission_;
    }

    friend std::ostream &operator<<(std::ostream &os, const FissionBank &bank);

private:
    const CoreMesh &mesh_;
    std::vector<Particle> sites_;
    real_t total_fission_;
};
} // namespace mc
} // namespace mocc
