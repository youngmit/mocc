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

#include <utility>
#include "util/pugifwd.hpp"
#include "core/core_mesh.hpp"
#include "core/solver.hpp"
#include "mc/fission_bank.hpp"
#include "mc/particle_pusher.hpp"

namespace mocc {
namespace mc {
/**
 * \brief Monte Carlo Eigenvalue solver
 *
 * This is essentially a driver for \ref mc::ParticlePusher, which performs
 * power iteration on fission banks. This class manages the particles in the
 * source/fission bank and maintains tallies for k-effective. Spatial tallies,
 * such as scalar flux and pin power are maintained within the \ref
 * mc::ParticlePusher. These tallies maintain batch statistics for each cycle,
 * which get reset at the end of the inactive cycles.
 */
class MonteCarloEigenvalueSolver : public Solver {
public:
    MonteCarloEigenvalueSolver(const pugi::xml_node &input,
                               const CoreMesh &mesh);

    /**
     * \brief Solve the eigenvalue problem
     *
          */
    void solve() override;

    /**
     * \brief Perform a single power iteration cycle
     */
    void step() override;

    void output(H5Node &node) const override;

private:
    // Data
    const CoreMesh &mesh_;
    const XSMesh xs_mesh_;
    ParticlePusher pusher_;
    int n_cycles_;
    int n_inactive_cycles_;
    int particles_per_cycle_;

    unsigned long seed_;

    RNG_LCG rng_;

    FissionBank source_bank_;

    // Cycle-by-cycle k history
    VecF k_history_tl_;
    VecF k_history_col_;
    VecF k_history_analog_;
    // Averaged k history
    VecF k_mean_history_tl_;
    VecF k_mean_history_col_;
    VecF k_mean_history_analog_;
    // K standard deviation history
    VecF k_stdev_history_tl_;
    VecF k_stdev_history_col_;
    VecF k_stdev_history_analog_;
    // Source bank Shannon entropy
    VecF h_history_;

    bool active_cycle_;

    std::pair<real_t, real_t> k_eff_;

    // Tally of results from the pusher_ tally, for computing batch statistics
    TallyScalar k_tally_tl_;
    TallyScalar k_tally_col_;
    TallyScalar k_tally_analog_;

    int cycle_;
    bool dump_sites_;
};
} // namespace mc
} // namespace mocc
