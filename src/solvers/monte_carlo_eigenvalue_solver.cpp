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
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/utils.hpp"
#include "mc/fission_bank.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
namespace mc {
MonteCarloEigenvalueSolver::MonteCarloEigenvalueSolver(
    const pugi::xml_node &input, const CoreMesh &mesh)
    : mesh_(mesh),
      xs_mesh_(mesh, MeshTreatment::TRUE),
      pusher_(mesh_, xs_mesh_),
      n_cycles_(input.attribute("cycles").as_int(-1)),
      n_inactive_cycles_(input.attribute("inactive_cycles").as_int(-1)),
      particles_per_cycle_(input.attribute("particles_per_cycle").as_int(-1)),
      seed_(input.attribute("seed").as_int(1)),
      rng_(seed_),
      source_bank_(input.child("fission_box"), particles_per_cycle_, mesh,
                   xs_mesh_, rng_),
      k_history_tl_(),
      k_history_col_(),
      k_history_analog_(),
      h_history_(),
      k_tally_tl_(),
      k_tally_col_(),
      k_tally_analog_(),
      cycle_(0),
      dump_sites_(false)
{
    // Check for valid input
    if (input.empty()) {
        throw EXCEPT("Input for Monte Carlo eigenvalue solver appears to "
                     "be empty.");
    }

    if (seed_ % 2 == 0) {
        throw EXCEPT("The RNG seed should be odd.");
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

    // Propagate the seed to the pusher
    pusher_.set_seed(seed_);

    return;
}
/**
 * The is pretty simple:
 *  - Loop over inactive cycles, calling step(),
 *  - clear the tally data on pusher_, then
 *  - loop over active cycles, calling step()
 */
void MonteCarloEigenvalueSolver::solve()
{
    cycle_ = -n_inactive_cycles_;
    LogScreen << "Performing inactive cycles:" << std::endl;
    LogScreen << std::setw(10) << "Cycle";
    LogScreen << std::setw(15) << "K-eff (TL)";
    LogScreen << std::setw(15) << "Mean (TL)";
    LogScreen << std::setw(15) << "Std. Dev. (TL)";
    LogScreen << std::setw(15) << "Mean (col)";
    LogScreen << std::setw(15) << "Mean (analog)";
    LogScreen << std::endl;
    active_cycle_ = false;
    for (int i = 0; i < n_inactive_cycles_; i++) {
        this->step();
    }

    // Reset the tallies following inactive cycles, here we want to reset ALL
    // tallies on the pusher_, scalar and spatial
    pusher_.reset_tallies(true);

    LogScreen << "Starting active cycles:" << std::endl;
    active_cycle_ = true;

    int n_active = n_cycles_ - n_inactive_cycles_ + 1;
    for (int i = 0; i < n_active; i++) {
        this->step();
    }

    return;
} // MonteCarloEigenvalueSolver::solve()

/**
 * Simulate all of the particles in the source bank using pusher_. After
 * simulating the batch of particles, extract eigenvalue estimates, and if in
 * active cycles, contribute to their tallies. At the end, swap storage between
 * the pusher_ fission bank and source_bank_, then resize source_bank_ to the
 * desired number of particles per cycle.
 */
void MonteCarloEigenvalueSolver::step()
{
    cycle_++;

    // Simulate all of the particles in the current fission bank
    pusher_.simulate(source_bank_, 1.0);

    // Log data
    k_eff_        = pusher_.k_tally_tl().get();
    auto k_tl     = pusher_.k_tally_tl().get();
    auto k_col    = pusher_.k_tally_col().get();
    auto k_analog = pusher_.k_tally_analog().get();

    k_history_tl_.push_back(k_tl.first);
    k_history_col_.push_back(k_col.first);
    k_history_analog_.push_back(k_analog.first);

    h_history_.push_back(source_bank_.shannon_entropy());

    if (active_cycle_) {
        k_tally_tl_.score(k_tl.first);
        k_tally_tl_.add_weight(1.0);
        k_tally_col_.score(k_col.first);
        k_tally_col_.add_weight(1.0);
        k_tally_analog_.score(k_analog.first);
        k_tally_analog_.add_weight(1.0);

        k_mean_history_tl_.push_back(k_tally_tl_.get().first);
        k_stdev_history_tl_.push_back(k_tally_tl_.get().second);
        k_mean_history_col_.push_back(k_tally_col_.get().first);
        k_stdev_history_col_.push_back(k_tally_col_.get().second);
        k_mean_history_analog_.push_back(k_tally_analog_.get().first);
        k_stdev_history_analog_.push_back(k_tally_analog_.get().second);
    }

    LogScreen << std::setw(10) << cycle_ << std::setw(15) << k_eff_.first;
    if (active_cycle_) {
        LogScreen << std::setw(15) << k_tally_tl_.get().first;
        LogScreen << std::setw(15) << k_tally_tl_.get().second;
        LogScreen << std::setw(15) << k_tally_col_.get().first;
        LogScreen << std::setw(15) << k_tally_analog_.get().first;
    }
    LogScreen << std::endl;

    // Grab the new fission sites from the pusher, and resize
    source_bank_.swap(pusher_.fission_bank());

    // Sort and re-index the source bank. This gives reproduceable IDs for
    // all particles, and therefore reproduceable parallel results. The
    // stable_sort is important.
    std::stable_sort(source_bank_.begin(), source_bank_.end());
    source_bank_.resize(particles_per_cycle_, rng_);
    unsigned i = 0;
    for (auto &p : source_bank_) {
        p.id = i++;
    }

    if (dump_sites_) {
        std::stringstream fname;
        fname << "sites_" << cycle_;
        std::ofstream f(fname.str());
        f << source_bank_;
    }

    // Reset the tallies on the Pusher. This should only reset the k-effective
    // tally, since all others are managed internally
    pusher_.reset_tallies();

    return;
} // MonteCarloEigenvalueSolver::step()

void MonteCarloEigenvalueSolver::output(H5Node &node) const
{
    auto dims = mesh_.dimensions();
    std::reverse(dims.begin(), dims.end());

    node.write("k_history_tl", k_history_tl_);
    node.write("k_history_col", k_history_col_);
    node.write("k_history_analog", k_history_analog_);

    node.write("h_history", h_history_);

    node.write("k_mean_history_tl", k_mean_history_tl_);
    node.write("k_mean_history_col", k_mean_history_col_);
    node.write("k_mean_history_analog", k_mean_history_analog_);

    node.write("k_stdev_history_tl", k_stdev_history_tl_);
    node.write("k_stdev_history_col", k_stdev_history_col_);
    node.write("k_stdev_history_analog", k_stdev_history_analog_);

    node.write("seed", seed_);

    pusher_.output(node);

    return;
}
} // namespace mc
} // namespace mocc
