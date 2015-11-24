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
        CMFD( const Mesh *mesh, SP_XSMeshHomogenized_t xsmesh );

        void solve( real_t &k, const ArrayB2 &flux );

        void project( ArrayX &flux );

        /**
         * Return a pointer to the coarse data. This is used to couple sweepers
         * and other objects that need access to the coarse data to the CMFD
         * solver.
         */
        CoarseData* get_data(){
            return &coarse_data_;
        }
    private:
        // Private methods
        void solve_1g( int group );

        // PRivate data
        const Mesh* mesh_;
        SP_XSMeshHomogenized_t xsmesh_;
        CoarseData coarse_data_;

        // Source vector
        Eigen::VectorXd b_;

        Source source_; 

        // One-group sparse matrix
        Eigen::SparseMatrix<real_t> m_;
    };
    typedef std::unique_ptr<CMFD> UP_CMFD_t;
}
