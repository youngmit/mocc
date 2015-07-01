#include "fixed_source_solver.hpp"

#include "transport_sweeper_factory.hpp"

namespace mocc{
FixedSourceSolver :: FixedSourceSolver( const pugi::xml_node &input ):
    sweeper_(std::move(TransportSweeperFactory(input)))
{

}

void FixedSourceSolver :: solve(){
}

}


