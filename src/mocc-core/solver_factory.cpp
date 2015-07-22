#include "solver_factory.hpp"

#include "error.hpp"

namespace mocc {
    SP_Solver_t SolverFactory( const pugi::xml_node &input,
        const CoreMesh &mesh ){

        if( input.empty() ) {
            Error("No input specified for the solver.");
        }
        return std::make_shared<EigenSolver>( input, mesh );
    }
}
