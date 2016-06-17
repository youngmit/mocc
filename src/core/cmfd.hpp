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

#include "coarse_data.hpp"
#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "mesh.hpp"
#include "source_isotropic.hpp"
#include "timers.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class CMFD {
    public:
        CMFD( const pugi::xml_node &input, const Mesh *mesh,
                SP_XSMeshHomogenized_t xsmesh );

        void solve( real_t &k );

        /**
         * \brief Return a pointer to the coarse data.
         *
         * This is used to couple sweepers and other objects that need access to
         * the coarse data to the CMFD solver.
         */
        CoarseData* get_data() {
            return &coarse_data_;
        }

        /**
         * \brief Return a reference to the \ref CoarseData object.
         */
        CoarseData& coarse_data() {
            return coarse_data_;
        }

        /**
         * \brief Return a reference to the MG flux.
         */
        const ArrayB2& flux() const {
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
        bool is_enabled() const {
            return is_enabled_;
        }

        /**
         * \brief Set the eigenvalue convergence tolerance
         */
        void set_k_tolerance( real_t tol ) {
            assert( tol > REAL_FUZZ );

            k_tol_ = tol;
            return;
        }

        /**
         * \brief Set the fission source convergence tolerance
         */
        void set_psi_tolerance( real_t tol ) {
            assert( tol > REAL_FUZZ );

            psi_tol_ = tol;
            return;
        }

    private:
        // Private methods
        void solve_1g( int group );
        void fission_source( real_t k );
        void print( int iter, real_t k, real_t k_err, real_t psi_err );

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

        const Mesh* mesh_;
        SP_XSMeshHomogenized_t xsmesh_;
        int n_cell_;
        CoarseData coarse_data_;
        bool is_enabled_;

        // Single-group fission source
        ArrayB1 fs_;
        ArrayB1 fs_old_;

        // Single-group flux guess
        VectorX x0_;

        SourceIsotropic source_;

        // Vector of one-group sparse matrix
        std::vector< Eigen::SparseMatrix<real_t> > m_;

        // Vector of BiCGSTAB objects.
        std::vector< Eigen::BiCGSTAB< Eigen::SparseMatrix<real_t> > > solvers_;

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
        int max_iter_;

        // Other options
        bool zero_fixup_;

    };
    typedef std::unique_ptr<CMFD> UP_CMFD_t;
}
