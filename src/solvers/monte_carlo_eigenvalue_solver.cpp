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

#include "monte_carlo_eigenvalue_solver.hpp"

#include <iomanip>

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"

#include "mc/fission_bank.hpp"

namespace {
const int WIDTH = 15;
}

namespace mocc {
namespace mc {
MonteCarloEigenvalueSolver::MonteCarloEigenvalueSolver(
    const pugi::xml_node &input, const CoreMesh &mesh)
    : mesh_(mesh),
      xs_mesh_(mesh),
      pusher_(mesh_, xs_mesh_),
      n_cycles_(input.attribute("cycles").as_int(-1)),
      n_inactive_cycles_(input.attribute("inactive_cycles").as_int(-1)),
      particles_per_cycle_(input.attribute("particles_per_cycle").as_int(-1)),
      source_bank_(input.child("fission_box"), particles_per_cycle_, mesh,
                   xs_mesh_),
      k_history_(),
      h_history_()
{
    // Check for valid input
    if (input.empty()) {
        throw EXCEPT("Input for Monte Carlo eigenvalue solver appears to "
                     "be empty.");
    }
    if (n_cycles_ < 0) {
        throw EXCEPT("Invalid number of cycle specified");
    }
    if (n_cycles_ == 0) {
        Warn("Zero cycles requested. You sure?");
    }

    if (n_inactive_cycles_ < 0) {
        throw EXCEPT("Invalid number of inactive cycles specified");
    }
    if (n_inactive_cycles_ == 0) {
        Warn("Zero inactive cycles requested. You sure?");
    }

    if (particles_per_cycle_ < 0) {
        throw EXCEPT("Invalid number of particles per cycle specified");
    }
    if (particles_per_cycle_ == 0) {
        Warn("Zero particles per cycle requested. You sure?");
    }

    return;
}

void MonteCarloEigenvalueSolver::solve()
{
    LogScreen << "Performing inactive cycles:" << std::endl;
    LogScreen << std::setw(15) << "Cycle K-eff";
    LogScreen << std::setw(15) << "Avg. K-eff";
    LogScreen << std::setw(15) << "Std. Dev.";
    LogScreen << std::endl;
    k_eff_        = {1.0, 0.0};
    active_cycle_ = false;
    for (int i = 0; i < n_inactive_cycles_; i++) {
        this->step();
    }

    LogScreen << "Starting active cycles:" << std::endl;
    active_cycle_ = true;

    int n_active = n_cycles_ - n_inactive_cycles_;
    for (int i = 0; i < n_active; i++) {
        this->step();
    }

    return;
} // MonteCarloEigenvalueSolver::solve()

void MonteCarloEigenvalueSolver::step()
{
    // Simulate all of the particles in the current fission bank
    pusher_.simulate(source_bank_, k_eff_.first);

    // Log data
    k_eff_ = pusher_.k_tally().get();
    k_history_.push_back(k_eff_.first);
    h_history_.push_back(source_bank_.shannon_entropy());
    if (active_cycle_) {
        k_tally_.score(k_eff_.first);
        k_tally_.add_weight(1.0);
        k_mean_history_.push_back(k_tally_.get().first);
        k_stdev_history_.push_back(k_tally_.get().second);
    }

    LogScreen << std::setw(15) << k_eff_.first;
    if (active_cycle_) {
        LogScreen << std::setw(15) << k_tally_.get().first;
        LogScreen << std::setw(15) << k_tally_.get().second;
    }
    LogScreen << std::endl;

    // Grab the new fission sites from the pusher, and resize
    source_bank_.swap(pusher_.fission_bank());
    source_bank_.resize(particles_per_cycle_);

    return;
} // MonteCarloEigenvalueSolver::step()

void MonteCarloEigenvalueSolver::output(H5Node &node) const
{
    node.write("k_history", k_history_);
    node.write("h_history", h_history_);
    node.write("k_mean_history", k_mean_history_);
    node.write("k_stdev_history", k_stdev_history_);
    node.create_group("flux");
    for(unsigned ig=0; ig<pusher_.flux_tallies().size(); ig++ ) {
        std::stringstream path;
        path << "flux/" << std::setfill('0') << std::setw(3) << ig;
        auto flux_result = pusher_.flux_tallies()[ig].get();
        VecF flux_val; flux_val.reserve(flux_result.size());
        VecF flux_stdev; flux_stdev.reserve(flux_result.size());
        for( const auto &v: flux_result ) {
            flux_val.push_back(v.first);
            flux_stdev.push_back(v.second);
        }
        node.write(path.str(), flux_val);
        path << "_stdev";
        node.write(path.str(), flux_stdev);
    }
    return;
}
} // namespace mc
} // namespace mocc
