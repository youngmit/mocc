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

#include <vector>

#include "core/core_mesh.hpp"
#include "core/output_interface.hpp"
#include "core/xs_mesh.hpp"

#include "fission_bank.hpp"
#include "particle.hpp"
#include "tally_scalar.hpp"
#include "tally_spatial.hpp"

namespace mocc {
namespace mc {
/**
 * This class manages the simulation of particle histories. Each call to
 * \ref simulate() will track the entire history of a particle until its
 * death and the death of all of its progeny (which would arise from
 * variance reduction techniques such as russian roulette/splitting).
 */
class ParticlePusher : public HasOutput {
public:
    ParticlePusher(const CoreMesh &mesh, const XSMesh &xs_mesh);

    /**
     * \brief Simulate a particle history
     *
     * \param p the \ref Particle to simulate
     * \param tally whether to treat the passed particle as a realization. This
     * should be true for history-based statistics, false for something else
     * like batch statistics.
     */
    void simulate(Particle p, bool tally=false);

    /**
     * \brief Simulate all particles in a \ref FissionBank
     *
     * \param bank the \ref FissionBank to generate/simulate particles from
     * \param k_eff the guess to use for k-effective to scale neutron
     * production in fission.
     */
    void simulate(const FissionBank &bank, real_t k_eff);

    /**
     * \brief Perform an interaction of a particle with its underlying
     * medium
     */
    void collide(Particle &p);

    /**
     * \brief Return a reference to the internal \ref FissionBank
     */
    FissionBank &fission_bank()
    {
        return fission_bank_;
    }

    /**
     * \brief Return a reference to the internal eigenvalue tally
     */
    const TallyScalar &k_tally() const
    {
        return k_tally_collision_;
    }

    /**
     * \brief Reset all tallies
     */
    void reset_tallies(bool clear_persistent = false)
    {
        k_tally_.reset();
        k_tally_collision_.reset();

        if (clear_persistent) {
            for (auto &flux_tally : scalar_flux_tally_) {
                flux_tally.reset();
            }

            for (auto &flux_tally : fine_flux_tally_) {
                flux_tally.reset();
            }

            for (auto &tally : fine_flux_col_tally_) {
                tally.reset();
            }
            pin_power_tally_.reset();
        }
        return;
    }

    void commit_tallies()
    {
        for (auto &t : scalar_flux_tally_) {
            t.commit_realization();
        }
        for (auto &t : fine_flux_tally_) {
            t.commit_realization();
        }
        for (auto &t : fine_flux_col_tally_) {
            t.commit_realization();
        }

        pin_power_tally_.commit_realization();

        return;
    }

    /**
     * \brief Assign a new seed to the RNG
     */
    void set_seed(unsigned long seed)
    {
        seed_ = seed;
    }

    const auto &flux_tallies() const
    {
        return scalar_flux_tally_;
    }

    const auto &fine_flux_tallies() const
    {
        return fine_flux_tally_;
    }

    void output(H5Node &node) const override;

private:
    const CoreMesh &mesh_;
    const XSMesh &xs_mesh_;

    int n_group_;

    // This fission bank stores new fission sites generated as the result of
    // simulating particles. This bank is cleared every time
    // simulate(FissionBank) is called
    FissionBank fission_bank_;

    // A map from mesh regions to XSMesh regions. This is useful since the
    // XSMesh goes from an XSMeshRegion to the mesh regions. It is more
    // necessary for MC to look up the cross sections for the current
    // region, somewhat at random.
    std::vector<int> xsmesh_regions_;

    // Do implicit capture?
    bool do_implicit_capture_;

    unsigned long seed_;

    // Eigenvalue tally
    TallyScalar k_tally_;
    TallyScalar k_tally_collision_;
    // Guess to use for scaling fission neutron production. Warning: Don't
    // try to use this as the actual system eigenvalue, since it is not tied
    // directly to a specific tally
    real_t k_eff_;

    // Scalar flux tallies
    std::vector<TallySpatial> scalar_flux_tally_;
    std::vector<TallySpatial> fine_flux_tally_;
    std::vector<TallySpatial> fine_flux_col_tally_;

    // Pin power tally
    TallySpatial pin_power_tally_;

    // Used to generate unique particle IDs
    unsigned id_offset_;

    unsigned n_cycles_;
    bool print_particles_;
};
} // namespace mc
} // namespace mocc
