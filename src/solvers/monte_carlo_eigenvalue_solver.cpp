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

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"

#include "mc/fission_bank.hpp"

namespace mocc {
MonteCarloEigenvalueSolver::MonteCarloEigenvalueSolver(
    const pugi::xml_node &input, const CoreMesh &mesh)
    : mesh_(mesh),
      xs_mesh_(mesh),
      pusher_(mesh_, xs_mesh_),
      n_cycles_(input.attribute("cycles").as_int(-1)),
      n_inactive_cycles_(input.attribute("inactive_cycles").as_int(-1)),
      particles_per_cycle_(input.attribute("particles_per_cycle").as_int(-1)),
      fission_bank_(input.child("fission_box"), particles_per_cycle_, mesh,
                    xs_mesh_) {
    // Check for valid input
    if (input.empty()) {
        throw EXCEPT(
            "Input for Monte Carlo eigenvalue solver appears to "
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

void MonteCarloEigenvalueSolver::solve() {
    k_eff_ = 0.0;
    k_variance_ = 0.0;

    LogScreen << "Performing inactive cycles:" << std::endl;
    active_cycle_ = false;
    for (int i = 0; i < n_inactive_cycles_; i++) {
        this->step();
    }

    LogScreen << "Starting active cycles:" << std::endl;
    active_cycle_ = true;
    for (int i = n_inactive_cycles_; i < n_cycles_; i++) {
        this->step();
    }

    return;
}  // MonteCarloEigenvalueSolver::solve()

void MonteCarloEigenvalueSolver::step() {
    // Simulate all of the particles in the current fission bank
//std::cout << fission_bank_;
std::cout << fission_bank_.size() << std::endl;
    pusher_.simulate(fission_bank_);


    real_t k_eff =
        pusher_.fission_bank().total_fission() / particles_per_cycle_;

    LogScreen << "K-effective: " << k_eff << std::endl;

    // Grab the new fission sites from the pusher, and resize
    fission_bank_.swap(pusher_.fission_bank());
std::cout << fission_bank_.size() << std::endl;
    fission_bank_.resize(particles_per_cycle_);

    return;
}  // MonteCarloEigenvalueSolver::step()

void MonteCarloEigenvalueSolver::output(H5Node &node) const { return; }

}  // namespace mocc
