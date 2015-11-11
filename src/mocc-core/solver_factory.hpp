#pragma once

#include "pugixml.hpp"

#include "mocc-core/eigen_solver.hpp"
#include "mocc-core/solver.hpp"

namespace mocc {
    /**
     * In a more complicated world, this would interrogate the input XML to
     * determine the type of highest-level \ref Solver to use, allocate and
     * construct that \ref Solver and return a shared pointer to it. In
     * practice, we are only making instances of \ref EigenSolver.
     */
    SP_Solver_t SolverFactory( const pugi::xml_node &input, 
        const CoreMesh &mesh );
}
