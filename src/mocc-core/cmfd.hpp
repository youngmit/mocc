#pragma once

#include <memory>

#include <Eigen/Sparse>

#include "coarse_data.hpp"
#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "mesh.hpp"
#include "source.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class CMFD {
    public:
        CMFD( const pugi::xml_node &input, const Mesh *mesh,
                SP_XSMeshHomogenized_t xsmesh );

        void solve( real_t &k, const ArrayB2 &flux );

        /**
         * Return a pointer to the coarse data. This is used to couple sweepers
         * and other objects that need access to the coarse data to the CMFD
         * solver.
         */
        CoarseData* get_data() {
            return &coarse_data_;
        }

        /**
         * Return a reference to the \ref CoarseData object.
         */
        CoarseData& coarse_data() {
            return coarse_data_;
        }

        /**
         * Return a referebce to the MG flux.
         */
        const ArrayB2& flux() const {
            return coarse_data_.flux;
        }

    private:
        // Private methods
        void solve_1g( int group );
        void fission_source( real_t k );

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
        const Mesh* mesh_;
        SP_XSMeshHomogenized_t xsmesh_;
        int n_cell_;
        CoarseData coarse_data_;

        // Single-group fission source
        ArrayF fs_;
        ArrayF fs_old_;

        Source source_;

        // Vector of one-group sparse matrix
        std::vector< Eigen::SparseMatrix<real_t> > m_;

        // Surface quantities. We need to keep these around to do the current
        // update without having to recalculate. Based on profiling, might be
        // nice to still get these on the fly to save on memory, but this is
        // fine for now.
        ArrayB2 d_hat_;
        ArrayB2 d_tilde_;
        ArrayB2 s_hat_;
        ArrayB2 s_tilde_;

        // Convergence options
        real_t k_tol_;
        real_t psi_tol_;
        int max_iter_;

    };
    typedef std::unique_ptr<CMFD> UP_CMFD_t;
}
