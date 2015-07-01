#pragma once

#include "pugixml.hpp"

#include "solver.hpp"
#include "eigen_solver.hpp"

namespace mocc {
    // In a more complicated world, this would interrogate the input XML to
    // determine the type of highest-level solver to use, allocate and construct
    // that solver and return a shared pointer to it. In practice, we are only
    // making eigenvalue solvers.
    SP_Solver_t SolverFactory( const pugi::xml_node &input );
}
