#include "solver_factory.hpp"

#include <string>

#include "pugixml.hpp"

#include "mocc-core/error.hpp"

namespace mocc {
    SP_Solver_t SolverFactory( const pugi::xml_node &input,
        const CoreMesh &mesh ){

        if( input.empty() ) {
            throw EXCEPT("No input specified for the solver.");
        }
        std::string type = input.attribute("type").value();
        if( type == "eigenvalue" ) {
            return std::make_shared<EigenSolver>( input, mesh );
        } else if( type == "fixed_source" ) {
            return std::make_shared<FixedSourceSolver>( input, mesh );
        } else {
            throw EXCEPT("Unrecognized solver type.");
        }
    }
}
