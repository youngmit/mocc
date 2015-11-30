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
         * Set up the linear systems for each group. This doesnt need to be done
         * for each iteration, nor in the case of non-zero currents/d-hat terms
         * can it, since the flux is allowed to change, which in turn will
         * affect the new D-hats. While this conusmes more memory to store the
         * systems for each group, it should be faster.
         */
        void setup_solve();
        real_t total_fission();

        // PRivate data
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

        // Convergence options
        real_t k_tol_;
        real_t psi_tol_;
        int max_iter_;

    };
    typedef std::unique_ptr<CMFD> UP_CMFD_t;
}
