#include "solver_factory.hpp"

namespace mocc {
    SP_Solver_t SolverFactory( pugi::xml_node &input ){
        return std::make_shared<EigenSolver>( input );
    }
}
