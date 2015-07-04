#include "fixed_source_solver.hpp"

#include "transport_sweeper_factory.hpp"

namespace mocc{
FixedSourceSolver :: FixedSourceSolver( const pugi::xml_node &input, 
        const CoreMesh &mesh )
{
    sweeper_ = UP_Sweeper_t((std::move(TransportSweeperFactory(input, mesh))));
}

void FixedSourceSolver :: solve(){
}

}


