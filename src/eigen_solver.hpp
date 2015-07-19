#pragma once

#include "pugixml.hpp"

#include "eigen_interface.hpp"
#include "solver.hpp"
#include "core_mesh.hpp"
#include "fixed_source_solver.hpp"

namespace mocc{

    class EigenSolver: public Solver{
    public:
        EigenSolver( const pugi::xml_node &input, const CoreMesh &mesh );
        void solve();
        void step();
    
    private:
        FixedSourceSolver fss_;
        
        // Fission source, and previous iterate
        MatrixX fission_source_;
        MatrixX fission_source_prev_;

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
    };
}
