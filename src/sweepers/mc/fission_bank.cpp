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

#include "fission_bank.hpp"

#include <iostream>
#include "pugixml.hpp"
#include "util/error.hpp"

namespace mocc {
namespace mc {
FissionBank::FissionBank(const CoreMesh &mesh) : mesh_(mesh)
{
    return;
}

/**
 * \brief Construct a FissionBank by uniformly sampling fission sites.
 *
 * \param input XML node containing bounds of a 3-D box within which to
 * sample initial fission sites
 * \param n the number of initial sites to sample
 * \param mesh the \ref CoreMesh to use for initial sampling
 * \param xs_mesh the \ref XSMesh to use for initial sampling
 *
 * This constructor initializes a \ref FissionBank using input specified in
 * an XML node.
 */
FissionBank::FissionBank(const pugi::xml_node &input, int n,
                         const CoreMesh &mesh, const XSMesh &xs_mesh,
                         RNG_LCG &rng)
    : mesh_(mesh), total_fission_(0.0)
{
    if (input.empty()) {
        throw EXCEPT("Empty input provided to FissionBank");
    }

    int ng = xs_mesh.n_group();

    // Make sure that all of the bounds are specified.
    if (input.attribute("x_min").empty() || input.attribute("x_max").empty() ||
        input.attribute("y_min").empty() || input.attribute("y_max").empty() ||
        input.attribute("z_min").empty() || input.attribute("z_max").empty()) {
        throw EXCEPT("Not all X, Y, Z bounds specified in fission_box");
    }

    real_t x_min = input.attribute("x_min").as_double(0.0);
    real_t x_max = input.attribute("x_max").as_double(0.0);
    real_t y_min = input.attribute("y_min").as_double(0.0);
    real_t y_max = input.attribute("y_max").as_double(0.0);
    real_t z_min = input.attribute("z_min").as_double(0.0);
    real_t z_max = input.attribute("z_max").as_double(0.0);

    // Make sure the bounds are valid
    if ((x_min >= x_max) || (y_min >= y_max) || (z_min >= z_max)) {
        throw EXCEPT("Invalid fission_box bounds specified.");
    }

    // See if we want to do a fissile region rejection (only accept fission
    // sites in fissile regions).
    bool fissile_rejection = input.attribute("fissile_rejection").as_bool(true);

    sites_.reserve(n);
    if (!fissile_rejection) {
        for (int i = 0; i < n; i++) {
            Point3 p(rng.random(x_min, x_max), rng.random(y_min, y_max),
                     rng.random(z_min, z_max));
            Direction dir(rng.random(TWOPI), rng.random(-HPI, HPI));
            int ig = rng.random_int(ng);
            sites_.emplace_back(p, dir, ig, i);
        }
    }
    else {
        Warn("Fissile region rejection is not supported yet.");
        for (int i = 0; i < n; i++) {
            Point3 p(rng.random(x_min, x_max), rng.random(y_min, y_max),
                     rng.random(z_min, z_max));
            Direction dir(rng.random(TWOPI), rng.random(-HPI, HPI));
            int ig = rng.random_int(ng);
            sites_.emplace_back(p, dir, ig, i);
        }
    }

    return;
}

real_t FissionBank::shannon_entropy() const
{
    real_t h = 0.0;

    VecF populations(mesh_.n_pin(), 0.0);

    for (const auto &p : sites_) {
        int icell = mesh_.coarse_cell_point(p.location_global);
        if (icell < 0 || icell > (int)mesh_.n_pin()) {
            Warn("ga");
        }
        populations[icell] += p.weight;
    }

    for (const auto &p : populations) {
        real_t pj = p / sites_.size();
        if (pj > 0.0) {
            h -= pj * std::log2(pj);
        }
    }

    return h;
}

void FissionBank::resize(unsigned int n, RNG_LCG &rng)
{
    assert(sites_.size() > 0);

    if (n > sites_.size()) {
        // Fission bank is too small. Randomly sample fission sites to
        // expand.
        // Make sure to only sample for the original sites. Probably not
        // absolutely necessary, but keeps the original sites equally
        // probable for the whole process.
        int n_orig = sites_.size();
        while (sites_.size() < n) {
            int i_rand = rng.random_int(n_orig);
            sites_.push_back(sites_[i_rand]);
        }
    }

    if (n < sites_.size()) {
        // Fission bank is too big. Randomly sample fission sites to remove.
        // We do a swap, then pop_back because A) the order doesnt matter,
        // and B) it avoids having to shuffle elements past the removed
        // element. I could probaly just truncate the vector, but for some
        // reason this feels more right
        int n_remove = sites_.size() - n;
        for (int i = 0; i < n_remove; i++) {
            int i_rand     = rng.random_int(sites_.size());
            sites_[i_rand] = sites_.back();
            sites_.pop_back();
        }
    }

    return;
}

void FissionBank::swap(FissionBank &other)
{
    sites_.swap(other.sites_);
    real_t tfis          = total_fission_;
    total_fission_       = other.total_fission_;
    other.total_fission_ = tfis;
}

std::ostream &operator<<(std::ostream &os, const FissionBank &bank)
{
    os << bank.sites_.size() << " fission sites:" << std::endl;
    for (const auto &p : bank.sites_) {
        os << p << std::endl;
    }
    return os;
}
} // namespace mc
} // namespace mocc
