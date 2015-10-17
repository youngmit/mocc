#pragma once

#include "pugixml.hpp"

#include "eigen_interface.hpp"
#include "solver.hpp"
#include "core_mesh.hpp"
#include "fixed_source_solver.hpp"
#include "transport_sweeper.hpp"
#include "h5file.hpp"
#include "cmfd.hpp"

namespace mocc{
    struct ConvergenceCriteria {
        ConvergenceCriteria( real_t k, real_t error_k, real_t error_psi ):
            k(k),
            error_k(error_k),
            error_psi(error_psi) { }
        real_t k;
        real_t error_k;
        real_t error_psi;
    };


    class EigenSolver: public Solver{
    public:
        EigenSolver( const pugi::xml_node &input, const CoreMesh &mesh );
        void solve();
        void step();
        
        const TransportSweeper* sweeper() const {
            return fss_.sweeper();
        }
        
        // Implement the output interface
        void output( H5::CommonFG *file ) const {
            VecF k;
            VecF error_k;
            VecF error_psi;

            for( auto &c: convergence_ ) {
                k.push_back(c.k);
                error_k.push_back(c.error_k);
                error_psi.push_back(c.error_psi);
            }

            VecI dims(1, convergence_.size());

            HDF::Write( file, "k", k, dims );
            HDF::Write( file, "error_k", error_k, dims );
            HDF::Write( file, "error_psi", error_psi, dims );

            fss_.output( file );
        }


    private:
        const static int out_w_ = 14;
        FixedSourceSolver fss_;
        
        // Fission source, and previous iterate
        ArrayF fission_source_;
        ArrayF fission_source_prev_;

        // Current guess for k
        real_t keff_;

        // Previous guess for k
        real_t keff_prev_;

        // Convergence criterion for the system eigenvalue
        real_t tolerance_k_;

        // Convergence criterion for the fission source distribution (L-2 norm)
        real_t tolerance_psi_;

        // Maximum allowable outer iterations
        unsigned int max_iterations_;

        // Print the current state of the eigenvalue solver
        void print( int iter, ConvergenceCriteria conv );

        // Vector of the convergence criteria. We will export these to the HDF5
        // file at the end of the run for posterity
        std::vector<ConvergenceCriteria> convergence_;

        // CMFD accelerator
        UP_CMFD_t cmfd_;
    };
}
