#include "eigen_solver.hpp"

namespace mocc{
    EigenSolver::EigenSolver( const pugi::xml_node &input, 
            const CoreMesh &mesh )
        : m_fss( input, mesh )
    {
    }

    void EigenSolver::solve() {
    }

    void EigenSolver::step() {

    }

};
