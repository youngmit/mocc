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

#include "cmfd.hpp"

#include <cmath>
#include <iomanip>
#include <vector>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/global_config.hpp"

typedef Eigen::Triplet<mocc::real_t> T;
typedef Eigen::SparseMatrix<mocc::real_t> M;

namespace mocc {
CMFD::CMFD(const pugi::xml_node &input, const Mesh *mesh,
           SP_XSMeshHomogenized_t xsmesh)
    : timer_(RootTimer.new_timer("CMFD", true)),
      timer_init_(timer_.new_timer("Initialization", true)),
      timer_setup_(timer_.new_timer("Setup Linear System")),
      timer_solve_(timer_.new_timer("Solve")),
      mesh_(mesh),
      xsmesh_(xsmesh),
      n_cell_(mesh->n_pin()),
      n_group_(xsmesh->n_group()),
      coarse_data_(*mesh, xsmesh->n_group()),
      is_enabled_(true),
      fs_(n_cell_),
      fs_old_(n_cell_),
      x_(n_cell_),
      source_(n_cell_, xsmesh_.get(), coarse_data_.flux),
      m_(n_group_, Eigen::SparseMatrix<real_t>(n_cell_, n_cell_)),
      solvers_(n_group_),
      d_hat_(mesh_->n_surf(), n_group_),
      d_tilde_(mesh_->n_surf(), n_group_),
      s_hat_(mesh_->n_surf(), n_group_),
      s_tilde_(mesh_->n_surf(), n_group_),
      n_solve_(0),
      k_tol_(1.0e-6),
      psi_tol_(1.0e-5),
      resid_reduction_(0.001),
      max_iter_(100),
      zero_fixup_(false)
{
    // Set up the structure of the matrix
    std::vector<T> structure;
    for (int i = 0; i < n_cell_; i++) {
        structure.push_back(T(i, i, 1.0));
        for (auto d : AllSurfaces) {
            int n = mesh_->coarse_neighbor(i, d);
            if (n >= 0) {
                // Defining both, though i could do this a little more
                // efficiently.
                structure.push_back(T(i, n, 1.0));
                structure.push_back(T(n, i, 1.0));
            }
        }
    }
    for (auto &m : m_) {
        m.setFromTriplets(structure.begin(), structure.end());
        m.makeCompressed();
    }

    // Parse options from the XML, if present
    if (!input.empty()) {
        // Eigenvalue tolerance
        if (!input.attribute("k_tol").empty()) {
            k_tol_ = input.attribute("k_tol").as_float(-1.0);
            if (k_tol_ <= 0.0) {
                throw EXCEPT("K tolerance is invalid.");
            }
        }

        // Fission source tolerance
        if (!input.attribute("psi_tol").empty()) {
            psi_tol_ = input.attribute("psi_tol").as_float(-1.0);
            if (psi_tol_ <= 0.0) {
                throw EXCEPT("Psi tolerance is invalid.");
            }
        }

        // Residual reduction
        if (!input.attribute("residual_reduction").empty()) {
            resid_reduction_ =
                input.attribute("residual_reduction").as_float(-1.0);
            if (resid_reduction_ <= 0.0) {
                throw EXCEPT("Residual reduction is invalid.");
            }
        }

        // Max iterations
        if (!input.attribute("max_iter").empty()) {
            max_iter_ = input.attribute("max_iter").as_int(-1);
            if (max_iter_ < 0) {
                throw EXCEPT("Max iterations invalid.");
            }
        }

        if (!input.attribute("negative_fixup").empty()) {
            zero_fixup_ = input.attribute("negative_fixup").as_bool(false);
        }

        // Enabled
        if (!input.attribute("enabled").empty()) {
            is_enabled_ = input.attribute("enabled").as_bool(true);
        }
    }

    timer_.toc();
    timer_init_.toc();
    return;
}

void CMFD::solve(real_t &k)
{
    // Make sure we have the requisite data
    if (!xsmesh_) {
        throw EXCEPT("No XS Mesh data! Need!");
    }
    timer_.tic();

    // Make sure no negative flux
    if (zero_fixup_) {
        for (int ig = 0; ig < (int)m_.size(); ig++) {
            ArrayB1 flux_1g = coarse_data_.flux(blitz::Range::all(), ig);
            // assert( flux_1g.isStorageContiguous() );
            for (int i = 0; i < (int)flux_1g.size(); i++) {
                flux_1g(i) = std::max(flux_1g(i), 0.0);
            }
        }
    }

    // Update homogenized cross sections
    xsmesh_->update();

    // Set up the linear systems
    this->setup_solve();

    timer_solve_.tic();

    real_t k_old = k;
    this->fission_source(k);
    real_t tfis = this->total_fission();

    // Calculate initial residual
    real_t r0 = this->residual();

    // Use the residual to set tolerance on the BiCGSTAB solvers
    for (auto &solver : solvers_) {
        solver.setTolerance(resid_reduction_ * r0);
    }

    auto flags = LogScreen.flags();
    LogScreen << "CMFD Converging to " << std::scientific << k_tol_ << " "
              << std::scientific << psi_tol_ << std::endl;
    LogScreen.flags(flags);

    int iter       = 0;
    real_t psi_err = 1.0;
    real_t ri      = 0.0; // Iteration residual
    while (true) {
        iter++;
        // Compute fission source
        fs_old_ = fs_;
        this->fission_source(k);
        real_t tfis_old = tfis;

        ri = 0.0;
        for (int group = 0; group < n_group_; group++) {
            source_.initialize_group(group);
            source_.fission(fs_, group);
            source_.in_scatter(group);
            source_.scale(mesh_->coarse_volume());

            this->solve_1g(group);
        }
        ri = std::sqrt(ri) / (n_cell_ * n_group_);

        tfis  = this->total_fission();
        k_old = k;
        k     = k * tfis / tfis_old;

        // Convergence check
        psi_err = 0.0;
        for (int i = 0; i < (int)fs_.size(); i++) {
            real_t e = fs_(i) - fs_old_(i);
            psi_err += e * e;
        }
        psi_err = std::sqrt(psi_err);

        if (((std::abs(k - k_old) < k_tol_) && (psi_err < psi_tol_) &&
             (ri / r0 < resid_reduction_)) ||
            (iter > max_iter_)) {
            break;
        }

        if ((iter % 10) == 0) {
            this->print(iter, k, std::abs(k - k_old), psi_err, ri / r0);
        }
    }
    this->print(iter, k, std::abs(k - k_old), psi_err, ri / r0);

    // Calculate the resultant currents and store back onto the coarse data
    this->store_currents();

    n_solve_++;

    timer_solve_.toc();
    timer_.toc();
    return;
} // solve()

real_t CMFD::solve_1g(int group)
{
    ArrayB1 flux_1g = coarse_data_.flux(blitz::Range::all(), group);

    for (int i = 0; i < n_cell_; i++) {
        x_[i] = flux_1g(i);
    }

    real_t resid = this->residual(group);

    x_ = solvers_[group].solveWithGuess(source_.get(), x_);

    // Store the result of the LS solution onto the CoarseData
    for (int i = 0; i < n_cell_; i++) {
        flux_1g(i) = x_[i];
    }

    return resid;
}

void CMFD::fission_source(real_t k)
{
    fs_ = 0.0;
    for (const auto &xsr : *xsmesh_) {
        for (const int i : xsr.reg()) {
            for (int ig = 0; ig < n_group_; ig++) {
                fs_(i) += xsr.xsmacnf(ig) * coarse_data_.flux(i, ig);
            }
        }
    }

    real_t r_keff = 1.0 / k;
    for (auto &v : fs_) {
        v *= r_keff;
    }
    return;
}

real_t CMFD::total_fission()
{
    real_t f = 0.0;
    for (const auto &xsr : *xsmesh_) {
        for (const int i : xsr.reg()) {
            for (int ig = 0; ig < n_group_; ig++) {
                f += xsr.xsmacnf(ig) * coarse_data_.flux(i, ig);
            }
        }
    }

    return f;
}

void CMFD::setup_solve()
{
    timer_setup_.tic();

    const Mesh::BCArray_t bc = mesh_->boundary_array();
    // Construct the system matrix
    size_t n_surf = mesh_->n_surf();
    int group     = 0;
    for (auto &m : m_) {
        // Diffusion coefficients
        VecF d_coeff(n_cell_);
        VecF xsrm(n_cell_);
        for (const auto &xsr : *xsmesh_) {
            real_t d  = 1.0 / (3.0 * xsr.xsmactr(group));
            real_t rm = xsr.xsmacrm(group);
            for (const int i : xsr.reg()) {
                d_coeff[i] = d;
                xsrm[i]    = rm;
            }
        }

        // Surface diffusivity (d_tilde) and non-linear correction
        // coefficient (d_hat)
        // There are lots of options to optimize this, mostly algebraic
        // simplifications, but this is very conformal to the canonical
        // formulations of CMFD found in the literature. If this starts
        // taking too much time, optimize.
        ArrayB1 d_tilde = d_tilde_(blitz::Range::all(), group);
        ArrayB1 d_hat   = d_hat_(blitz::Range::all(), group);
        ArrayB1 s_tilde = s_tilde_(blitz::Range::all(), group);
        ArrayB1 s_hat   = s_hat_(blitz::Range::all(), group);

        // Loop over the surfaces in the mesh, and calculate the inter-cell
        // coupling coefficients
        for (int is = 0; is < (int)n_surf; is++) {
            auto cells  = mesh_->coarse_neigh_cells(is);
            Normal norm = mesh_->surface_normal(is);

            real_t diffusivity_1 = 0.0;
            real_t diffusivity_2 = 0.0;
            if (cells.first > -1) {
                diffusivity_1 = d_coeff[cells.first] /
                                mesh_->cell_thickness(cells.first, norm);
            } else {
                switch (bc[(int)(norm)][0]) {
                case Boundary::REFLECT:
                    diffusivity_1 = 0.0 / 2.0;
                    break;
                case Boundary::VACUUM:
                    diffusivity_1 = 0.5 / 2.0;
                    break;
                default:
                    throw EXCEPT("Unsupported boundary type");
                }
            }

            if (cells.second > -1) {
                diffusivity_2 = d_coeff[cells.second] /
                                mesh_->cell_thickness(cells.second, norm);
            } else {
                switch (bc[(int)(norm)][1]) {
                case Boundary::REFLECT:
                    diffusivity_2 = 0.0 / 2.0;
                    break;
                case Boundary::VACUUM:
                    diffusivity_2 = 0.5 / 2.0;
                    break;
                default:
                    throw EXCEPT("Unsupported boundary type");
                }
            }

            d_tilde(is) = 2.0 * diffusivity_1 * diffusivity_2 /
                          (diffusivity_1 + diffusivity_2);

            // S-tilde is a mess. Since surface flux is calculated as
            // phi = s_tilde*flux_left + (1-s_tilde)*flux_right, there is an
            // inherent binding to a cell, as well as a surface. We assume
            // the convention that if possible the bound cell is the one to
            // the "left" of the surface. When such a cell is not present
            // (domain boundary), the cell is to the "right"
            real_t stil = (diffusivity_1 > 0.0)
                              ? diffusivity_1 / (diffusivity_1 + diffusivity_2)
                              : diffusivity_2 / (diffusivity_1 + diffusivity_2);

            s_tilde(is) = stil;

            // If we have currents defined from a transport sweeper or the
            // like, calculate D-hat coefficients
            bool have_data = norm == Normal::Z_NORM
                                 ? coarse_data_.has_axial_data()
                                 : coarse_data_.has_radial_data();
            if (have_data) {
                real_t j        = coarse_data_.current(is, group);
                real_t sfc_flux = coarse_data_.surface_flux(is, group);
                real_t flux_l   = cells.first >= 0
                                    ? coarse_data_.flux(cells.first, group)
                                    : 0.0;
                real_t flux_r = cells.second >= 0
                                    ? coarse_data_.flux(cells.second, group)
                                    : 0.0;
                d_hat(is) =
                    (j + d_tilde(is) * (flux_r - flux_l)) / (flux_l + flux_r);
                s_hat(is) = (cells.first >= 0)
                                ? (sfc_flux - s_tilde(is) * flux_l -
                                   (1.0 - s_tilde(is)) * flux_r) /
                                      (flux_l + flux_r)
                                : (sfc_flux - s_tilde(is) * flux_r) / (flux_r);
            } else {
                d_hat(is) = 0.0;
                s_hat(is) = 0.0;
            }
        } // surfaces

        // put values into the matrix. Optimal access patterns in sparse
        // matrix representations are not obvious, so the best way is to
        // iterate through the matrix linearly and act according to the
        // indices (i.e.  row/col) that we get for each coefficient.
        for (int k = 0; k < m.outerSize(); k++) {
            for (M::InnerIterator it(m, k); it; ++it) {
                auto i = it.row();
                auto j = it.col();
                if (i == j) {
                    // Diagonal element
                    real_t v = mesh_->coarse_volume(i) * xsrm[i];
                    for (auto is : AllSurfaces) {
                        size_t surf = mesh_->coarse_surf(i, is);
                        real_t a    = mesh_->coarse_area(i, is);
                        // Switch sign of D-hat if necessary
                        real_t d_hat_ij = d_hat(surf);
                        if ((is == Surface::WEST) || (is == Surface::SOUTH) ||
                            (is == Surface::BOTTOM)) {
                            d_hat_ij = -d_hat_ij;
                        }

                        v += a * (d_tilde(surf) + d_hat_ij);
                    }
                    it.valueRef() = v;
                } else {
                    // off-diagonal element
                    auto pair       = mesh_->coarse_interface(i, j);
                    real_t a        = mesh_->coarse_area(i, pair.second);
                    size_t surf     = pair.first;
                    real_t d_hat_ij = d_hat(surf);
                    // Switch sign of D-hat if necessary
                    if ((pair.second == Surface::WEST) ||
                        (pair.second == Surface::SOUTH) ||
                        (pair.second == Surface::BOTTOM)) {
                        d_hat_ij = -d_hat_ij;
                    }

                    it.valueRef() = a * (d_hat_ij - d_tilde(surf));
                }
            }
        } // matrix element loop

        solvers_[group].compute(m);
        solvers_[group].setMaxIterations(1500);

        group++;
    } // group loop
    timer_setup_.toc();
    return;
} // setup_solve

void CMFD::store_currents()
{
    int n_group = xsmesh_->n_group();
    int n_surf  = mesh_->n_surf();
    for (int ig = 0; ig < n_group; ig++) {
        auto all             = blitz::Range::all();
        auto current_1g      = coarse_data_.current(blitz::Range::all(), ig);
        auto surface_flux_1g = coarse_data_.surface_flux(all, ig);
        auto partial_1g      = coarse_data_.partial_current(all, ig);
        auto partial_old_1g  = coarse_data_.partial_current_old(all, ig);

        partial_old_1g = partial_1g;
        coarse_data_.set_has_old_partial(n_solve_ > 0);

        for (int is = 0; is < n_surf; is++) {
            auto cells = mesh_->coarse_neigh_cells(is);
            real_t flux_r =
                cells.second >= 0 ? coarse_data_.flux(cells.second, ig) : 0.0;
            real_t flux_l =
                cells.first >= 0 ? coarse_data_.flux(cells.first, ig) : 0.0;

            real_t d_hat   = d_hat_(is, ig);
            real_t d_tilde = d_tilde_(is, ig);
            real_t current =
                -d_tilde * (flux_r - flux_l) + d_hat * (flux_r + flux_l);

            current_1g(is) = current;

            real_t s_hat   = s_hat_(is, ig);
            real_t s_tilde = s_tilde_(is, ig);
            real_t surface_flux =
                (cells.first >= 0)
                    ? s_tilde * flux_l + (1.0 - s_tilde) * flux_r +
                          s_hat * (flux_l + flux_r)
                    : s_tilde * flux_r + s_hat * (flux_l + flux_r);

            surface_flux_1g(is) = surface_flux;

            partial_1g(is) = {{0.25 * surface_flux + 0.5 * current,
                               0.25 * surface_flux - 0.5 * current}};
        } // surfaces
    }     // groups
    return;
}

real_t CMFD::residual()
{
    real_t norm = 0.0;

    for (int group = 0; group < n_group_; group++) {
        ArrayB1 flux_1g = coarse_data_.flux(blitz::Range::all(), group);
        source_.initialize_group(group);
        source_.fission(fs_, group);
        source_.in_scatter(group);
        source_.scale(mesh_->coarse_volume());

        for (int icell = 0; icell < n_cell_; ++icell) {
            x_[icell] = flux_1g(icell);
        }

        norm += this->residual(group);
    }

    // int group = 0;
    // for (const auto &m : m_) {
    //    ArrayB1 flux_1g = coarse_data_.flux(blitz::Range::all(), group);
    //    VectorX flux_vec(n_cell_);
    //    // memcpy would be nice here...
    //    for (int icell = 0; icell < n_cell_; ++icell) {
    //        flux_vec[icell] = flux_1g(icell);
    //    }

    //    source_.initialize_group(group);
    //    source_.fission(fs_, group);
    //    source_.in_scatter(group);
    //    source_.scale(mesh_->coarse_volume());

    //    VectorX residual = m * flux_vec - source_.get();
    //    norm += residual.squaredNorm();

    //    group++;
    //}

    return std::sqrt(norm) / (n_group_ * n_cell_);
}

real_t CMFD::residual(int group) const
{
    VectorX resid = m_[group] * x_ - source_.get();

    return resid.squaredNorm();
}

void CMFD::print(int iter, real_t k, real_t k_err, real_t psi_err,
                 real_t resid_ratio)
{
    auto flags = LogScreen.flags();
    LogScreen << "       " << std::setprecision(5) << std::setw(6) << std::fixed
              << RootTimer.time() << " " << iter << " " << std::setprecision(10)
              << k << " " << std::scientific << k_err << " " << std::scientific
              << psi_err << " " << std::scientific << resid_ratio << std::endl;
    LogScreen.flags(flags);
    return;
}
}
