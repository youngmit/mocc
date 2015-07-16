#include "eigen_solver.hpp"

namespace mocc{
    EigenSolver::EigenSolver( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        fss_( input, mesh ),
        fission_source_( fss_.n_reg(), 1 )
    {
        
    }

    void EigenSolver::solve() {
    }

    void EigenSolver::step() {

    }

};
