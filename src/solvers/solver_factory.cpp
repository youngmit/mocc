#include "solver_factory.hpp"

#include <string>

#include "pugixml.hpp"

#include "core/error.hpp"
#include "core/files.hpp"

#include "eigen_solver.hpp"
#include "fixed_source_solver.hpp"

namespace mocc {
    SP_Solver_t SolverFactory( const pugi::xml_node &input,
        const CoreMesh &mesh ){

        LogFile << "Initializing solver..." << std::endl;

        SP_Solver_t solver;

        if( input.empty() ) {
            throw EXCEPT("No input specified for the solver.");
        }
        std::string type = input.attribute("type").value();
        if( type == "eigenvalue" ) {
            solver = std::make_shared<EigenSolver>( input, mesh );
        } else if( type == "fixed_source" ) {
            solver = std::make_shared<FixedSourceSolver>( input, mesh );
        } else {
            throw EXCEPT("Unrecognized solver type.");
        }

        LogFile << "Done initializing solver." << std::endl;

        return solver;
    }
}
