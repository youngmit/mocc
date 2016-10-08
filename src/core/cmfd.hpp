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

#include <memory>

#include <Eigen/Sparse>

#include "util/global_config.hpp"
#include "util/timers.hpp"
#include "coarse_data.hpp"
#include "eigen_interface.hpp"
#include "mesh.hpp"
#include "source_isotropic.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
class CMFD : public HasOutput {
public:
    CMFD(const pugi::xml_node &input, const CoreMesh *mesh,
         SP_XSMeshHomogenized_t xsmesh);

    /**
     * \brief Solve the CMFD system
     *
     * \param [in,out] k the initial guess to use for the system eigenvalue.
     * Updated byt the solve.
     */
    void solve(real_t &k);

    /**
     * \brief Return a pointer to the coarse data.
     *
     * This is used to couple sweepers and other objects that need access to
     * the coarse data to the CMFD solver.
     */
    CoarseData *get_data()
    {
        return &coarse_data_;
    }

    /**
     * \brief Return a reference to the \ref CoarseData object.
     */
    CoarseData &coarse_data()
    {
        return coarse_data_;
    }

    /**
     * \brief Return a reference to the MG flux.
     */
    const ArrayB2 &flux() const
    {
        return coarse_data_.flux;
    }

    /**
     * \brief Return whether the CMFD object is meant to actually do a
     * solve.
     *
     * It is possible to need a CMFD object to be constructed to use its
     * \ref CoarseData and other functionality, but not want to perform
     * actual CMFD solves.
     */
    bool is_enabled() const
    {
        return is_enabled_;
    }

    /**
     * \brief Set the eigenvalue convergence tolerance
     */
    void set_k_tolerance(real_t tol)
    {
        assert(tol > REAL_FUZZ);

        k_tol_ = tol;
        return;
    }

    /**
     * \brief Set the fission source convergence tolerance
     */
    void set_psi_tolerance(real_t tol)
    {
        assert(tol > REAL_FUZZ);

        psi_tol_ = tol;
        return;
    }

    void output(H5Node &node) const;

private:
    // Private methods
    /**
     * \brief Compute and return the L-2 norm of the residual of the CMFD
     * system, scaled by the size of the system.
     *
     * This calculates the residual of the CMFD system, considering all groups.
     * Since this implementation of the method must calculate the source for all
     * groups and needs to do repeated copies of the flux, it is somewhat more
     * expensive than the single-group implementation, which assumes that the
     * source is already calculated. As a result this is useful for calculating
     * the initial residual at the beginning of the CMFD solve, while the other
     * is more ideal for calculating the residual for each iteration.
     *
     * \pre The group-independent fission source has been calculated and stored
     * in \c fs_.
     *
     * \sa CMFD::residual(int)
     */
    real_t residual();

    /**
     * \brief Compute and return the squared L-2 norm of the residual for a
     * single group of the CMFD system
     *
     * \param group the group for which to calculate residual
     *
     * This is similar to the all-group method, but it only computes the
     * contribution to the residual for the requested group. It assumes that the
     * source (right-hand side) has already been computed and is therefore
     * cheaper to use for assessing convergence.
     *
     * Unlike \ref CMFD::residual(), this returns an un-scaled, squared L-2
     * norm. This should lead to the following code producing the same result as
     * the group-independent \ref CMFD::residual() method:
     * \code
norm = 0.0;
for(int ig=0; ig<n_group_; ig++) {
    // Satisfy preconditions
    norm += this->residual(group);
}
norm = std::sqrt(norm)/(n_cell_*n_group_);
     * \endcode
     *
     * \pre The group-independent fission source has been calculated and stored
     * in \c fs_.
     *
     * \pre The source (fission and inscattering) must be calculated for the
     * current group.
     *
     * \pre The current-group flux is already stored in the \c x_ vector
     *
     * \sa CMFD::residual()
     */
    real_t residual(int group) const;
    real_t solve_1g(int group);
    void fission_source(real_t k);
    void print(int iter, real_t k, real_t k_err, real_t psi_err,
               real_t resid_ratio);

    /**
     * \brief Calculate CMFD-derived currents after a CMFD solve and store
     * on the \ref CoarseData.
     *
     * After performing a CMFD solve, it is useful to have access to the
     * state of the current, as predicted by the CMFD system given the
     * current state of the D-hats. This allows us to calculate, for
     * instance, transverse leakage for an MoC sweeper, which otherwise
     * wouldn't know what the current looks like in the axial dimension.
     */
    void store_currents();

    /**
     * Set up the linear systems for each group. This doesnt need to be done
     * for each iteration, nor in the case of non-zero currents/d-hat terms
     * can it, since the flux is allowed to change, which in turn will
     * affect the new D-hats. While this conusmes more memory to store the
     * systems for each group, it should be faster.
     */
    void setup_solve();
    real_t total_fission();

    // Private data
    Timer &timer_;
    Timer &timer_init_;
    Timer &timer_setup_;
    Timer &timer_solve_;

    Mesh mesh_;
    const Mesh *fine_mesh_;
    XSMeshHomogenized xsmesh_;
    int n_cell_;
    int n_surf_;
    int n_group_;
    CoarseData coarse_data_;
    bool is_enabled_;

    // Single-group fission source
    ArrayB1 fs_;
    ArrayB1 fs_old_;

    // Single-group flux result from LS solve, also used for initial guess
    VectorX x_;

    SourceIsotropic source_;

    // Vector of one-group sparse matrix
    std::vector<Eigen::SparseMatrix<real_t>> m_;

    // Vector of BiCGSTAB objects.
    std::vector<Eigen::BiCGSTAB<Eigen::SparseMatrix<real_t>>> solvers_;

    // Surface quantities. We need to keep these around to do the current
    // update without having to recalculate. Based on profiling, might be
    // nice to still get these on the fly to save on memory, but this is
    // fine for now.
    ArrayB2 d_hat_;
    ArrayB2 d_tilde_;
    ArrayB2 s_hat_;
    ArrayB2 s_tilde_;

    // Number of times solve() has been called
    int n_solve_;

    // Convergence options
    real_t k_tol_;
    real_t psi_tol_;
    real_t resid_reduction_;
    int max_iter_;

    // Other options
    bool zero_fixup_;
    bool dump_current_;
};
typedef std::unique_ptr<CMFD> UP_CMFD_t;
}
