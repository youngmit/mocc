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

#include "eigen_solver.hpp"

#include <iomanip>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/string_utils.hpp"
#include "util/utils.hpp"
#include "util/validate_input.hpp"
#include "core/globals.hpp"

const static int out_w = 14;

namespace {
const std::vector<std::string> recognized_attributes = {
    "type", "cmfd", "k_tol", "psi_tol", "max_iter", "min_iter"};
}

namespace mocc {
EigenSolver::EigenSolver(const pugi::xml_node &input, const CoreMesh &mesh)
    : fss_(input, mesh),
      fission_source_(fss_.sweeper()->n_reg_fission()),
      fission_source_prev_(fss_.sweeper()->n_reg_fission()),
      min_iterations_(0)
{
    LogFile << "Initializing Eigenvalue solver..." << std::endl;

    validate_input(input, recognized_attributes);

    if (input.empty()) {
        throw EXCEPT("No input specified for the eigenvalue solver.");
    }

    // grab the convergence constraints from the XML
    int in_int      = 0;
    real_t in_float = 0.0;

    // K tolerance
    in_float = input.attribute("k_tol").as_float(-1.0);
    if (in_float <= 0.0) {
        throw EXCEPT("Invalid k tolerance.");
    }
    tolerance_k_ = in_float;

    // Psi tolerance
    in_float = input.attribute("psi_tol").as_float(-1.0);
    if (in_float <= 0.0) {
        throw EXCEPT("Invalid psi tolerance.");
    }
    tolerance_psi_ = in_float;

    // Max iterations
    in_int = input.attribute("max_iter").as_int(-1);
    if (in_int < 0) {
        throw EXCEPT("Invalid number of maximum iterations.");
    }
    max_iterations_ = in_int;

    // Min iterations
    if (!input.attribute("min_iter").empty()) {
        in_int = input.attribute("min_iter").as_int(-1);
        if ((in_int < 0) || (in_int > (int)max_iterations_)) {
            throw EXCEPT("Invalid number of minimum iterations.");
        }
        min_iterations_ = in_int;
    }

    // Read in dump iterations if present
    if (!input.child("dump_iterations").empty()) {
        dump_iterations_ =
            explode_string<int>(input.child("dump_iterations").child_value());
        std::sort(dump_iterations_.begin(), dump_iterations_.end());
        LogFile << "Dumping data for eigenvalue iterations:" << std::endl;
        for (const auto i : dump_iterations_) {
            LogFile << i << " ";
        }
        LogFile << std::endl;
    }

    // Count the number of fissile mesh regions
    n_fissile_regions_ = 0;
    for (const auto &xsr : fss_.sweeper()->xs_mesh()) {
        bool has_fission = false;
        for (int ig = 0; ig < xsr.n_group(); ig++) {
            if (xsr.xsmacnf(ig) > 0.0) {
                has_fission = true;
                break;
            }
        }
        if (has_fission) {
            n_fissile_regions_ += xsr.reg().size();
        }
    }

    // CMFD acceleration
    bool do_cmfd = input.attribute("cmfd").as_bool(false);
    if (do_cmfd) {
        // construct the CMFD solver using the mesh from the transport
        // sweeper
        cmfd_.reset(new CMFD(input.child("cmfd"), &mesh,
                             fss_.sweeper()->get_homogenized_xsmesh()));
        // Associate the sweeper with the coarse data from the CMFD solver
        CoarseData *const cd = cmfd_->get_data();
        fss_.sweeper()->set_coarse_data(cd);
    }

    LogFile << "Done initializing Eigenvalue solver." << std::endl;

    return;
}

// Perform a full-blown eigenvalue solve. Start with a guess for the
// fission source (flat), start doing power iteration. Once we have that
// working, we will start factoring in CMFD and other fancy tricks
void EigenSolver::solve()
{

    LogScreen << "Converging to: \n"
                 "\t Eigenvalue: "
              << tolerance_k_ << "\n"
              << "\t Fission Source (L-2 norm): " << tolerance_psi_ << "\n"
              << "\t Min/Max Iterations: " << min_iterations_ << " / "
              << max_iterations_ << "\n\n";

    keff_      = 1.0;
    keff_prev_ = 1.0;

    // initialize the fixed source solver and calculation the initial
    // fission source
    fss_.initialize();

    // Hand a reference to the fission source to the fixed source solver
    fss_.set_fission_source(&fission_source_);

    error_k_   = tolerance_k_;   // K residual
    error_psi_ = tolerance_psi_; // L-2 norm of the fission source residual

    fss_.sweeper()->calc_fission_source(keff_, fission_source_);

    LogScreen << std::setw(out_w) << "Time" << std::setw(out_w) << "Iter."
              << std::setw(out_w) << "k" << std::setw(out_w) << "k error"
              << std::setw(out_w) << "psi error" << std::endl;

    auto dump_it       = dump_iterations_.begin();
    unsigned next_dump = std::numeric_limits<unsigned>::max();
    if (dump_it != dump_iterations_.end()) {
        next_dump = *dump_it;
    }

    for (size_t n_iterations = 0; n_iterations < max_iterations_;
         n_iterations++) {
        this->step();

        // Check for convergence
        error_k_ = fabs(keff_ - keff_prev_);

        const auto &vol = fss_.sweeper()->volumes();
        //assert(vol.size() == fission_source_.size());
        Normalize(fission_source_.begin(), fission_source_.end());
        Normalize(fission_source_prev_.begin(),
                        fission_source_prev_.end());

        real_t efis = 0.0;
        for (int i = 0; i < (int)fss_.sweeper()->n_reg(); i++) {
            real_t e = (fission_source_(i) - fission_source_prev_(i));
            efis += e * e;
        }
        error_psi_ = std::sqrt(efis / n_fissile_regions_);

        convergence_.push_back(
            ConvergenceCriteria(keff_, error_k_, error_psi_));

        iteration_times_.push_back(RootTimer.time());

        this->print(n_iterations + 1, convergence_.back());

        if (n_iterations + 1 == next_dump) {
            std::stringstream fname;
            fname << global::case_name << "_iter_" << n_iterations + 1 << ".h5";
            H5Node h5f(fname.str(), H5Access::WRITE);
            this->output(h5f);

            next_dump = std::numeric_limits<int>::max();
            if (++dump_it != dump_iterations_.end()) {
                next_dump = *dump_it;
            }
        }

        // Check for NaN
        if (keff_ != keff_) {
            throw EXCEPT("Eigenvalue is not a number. Giving up.");
        }

        LogFile << (error_k_ < tolerance_k_) << " "
                << (error_psi_ < tolerance_psi_) << " "
                << (n_iterations >= min_iterations_) << std::endl;
        if ((error_k_ < tolerance_k_) && (error_psi_ < tolerance_psi_) &&
            (n_iterations >= min_iterations_)) {
            LogScreen << "Convergence criteria satisfied!" << std::endl;
            break;
        }

        if (n_iterations == (max_iterations_ - 1)) {
            LogScreen << "Maximum number of iterations reached!" << std::endl;
        }
    }
} // solve()

void EigenSolver::step()
{

    if (cmfd_ && cmfd_->is_enabled()) {
        this->do_cmfd();
    }

    // Perform a group sweep with the FSS
    // We are recalculating the fission source perhaps more than necessary.
    // Make if we want to remove redundant calculations, make sure that it will
    // jive with the normalization that is being done for the convergence
    // criterion.
    fss_.sweeper()->calc_fission_source(keff_, fission_source_);
    fission_source_prev_ = fission_source_;
    fss_.step();

    // Get the total fission sources
    real_t tfis1 = fss_.sweeper()->total_fission(false);
    real_t tfis2 = fss_.sweeper()->total_fission(true);

    // update estimate for k
    keff_prev_ = keff_;
    keff_      = keff_ * tfis1 / tfis2;

    // update the fission source
    fss_.sweeper()->calc_fission_source(keff_, fission_source_);

    return;
}

void EigenSolver::print(int iter, ConvergenceCriteria conv)
{
    LogScreen << std::setw(out_w) << std::fixed << std::setprecision(5)
              << RootTimer.time() << std::setw(out_w) << iter << conv
              << std::endl;
    return;
}

void EigenSolver::do_cmfd()
{
    assert(cmfd_);
    assert(cmfd_->is_enabled());
    // push homogenized flux onto the coarse mesh, solve, and pull it
    // back.
    cmfd_->coarse_data().flux =
        fss_.sweeper()->get_pin_flux(MeshTreatment::PIN_PLANE);

    // Set the convergence criteria for this solve, there are a few ways we
    // can do this. This needs some work. Should probably be converging the
    // CMFD-transport sweeper residual, rather than the actual CMFD solution.
    CMFDConvergence conv = CMFDConvergence::FIXED;
    switch (conv) {
    case CMFDConvergence::FIXED:
        break;
    case CMFDConvergence::FLOAT:
        real_t k_tol = std::max(error_k_ / 1000.0, tolerance_k_ / 10.0);
        cmfd_->set_k_tolerance(k_tol);

        real_t psi_tol = std::max(error_psi_ / 1000.0, tolerance_psi_ / 10.0);
        cmfd_->set_psi_tolerance(psi_tol);
        break;
    }
    cmfd_->solve(keff_);
    fss_.sweeper()->set_pin_flux(cmfd_->flux(), MeshTreatment::PIN_PLANE);
    //fss_.sweeper()->update_incoming_flux();
    return;
}

void EigenSolver::output(H5Node &file) const
{
    VecF k;
    VecF error_k;
    VecF error_psi;

    for (auto &c : convergence_) {
        k.push_back(c.k);
        error_k.push_back(c.error_k);
        error_psi.push_back(c.error_psi);
    }

    VecI dims(1, convergence_.size());

    {
        auto g = file.create_group("convergence");

        g.write("k", k, dims);
        g.write("error_k", error_k, dims);
        g.write("error_psi", error_psi, dims);
        g.write("iteration_time", iteration_times_);
        g.write("abscissae", iteration_times_);
    }

    fss_.output(file);
    if (cmfd_) {
        cmfd_->output(file);
    }
    return;
}

std::ostream &operator<<(std::ostream &os, ConvergenceCriteria conv)
{
    std::ios::fmtflags flags = os.flags();

    os << std::setw(out_w) << std::fixed << std::setprecision(10) << conv.k
       << std::setw(out_w) << std::scientific << std::setprecision(6)
       << conv.error_k << std::setw(out_w)
       << std::setiosflags(std::ios::scientific) << std::setprecision(6)
       << conv.error_psi;

    os.flags(flags);
    return os;
}
};
