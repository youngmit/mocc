#pragma once

#include "pugixml.hpp"

#include "eigen_interface.hpp"
#include "solver.hpp"
#include "core_mesh.hpp"
#include "fixed_source_solver.hpp"
#include "transport_sweeper.hpp"

namespace mocc{

    class EigenSolver: public Solver{
    public:
        EigenSolver( const pugi::xml_node &input, const CoreMesh &mesh );
        void solve();
        void step();
        
        const TransportSweeper* sweeper() const {
            return fss_.sweeper();
        }
    
    private:
        const static int out_w_ = 14;
        FixedSourceSolver fss_;
        
        // Fission source, and previous iterate
        ArrayX fission_source_;
        ArrayX fission_source_prev_;

        // Current guess for k
        float_t keff_;

        // Previous guess for k
        float_t keff_prev_;

        // Convergence criterion for the system eigenvalue
        float_t tolerance_k_;

        // Convergence criterion for the fission source distribution (L-2 norm)
        float_t tolerance_psi_;

        // Maximum allowable outer iterations
        unsigned int max_iterations_;

        // Print the current state of the eigenvalue solver
        void print( int iter, float_t error_k, float_t error_psi );
    };
}
