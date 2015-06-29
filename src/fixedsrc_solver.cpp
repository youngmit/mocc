#include "fixedsrc_solver.hpp"

#include "transport_sweeper_factory.hpp"

namespace mocc{
FixedSourceSolver :: FixedSourceSolver(pugi::xml_node &input):
    sweeper_(TransportSweeperFactory(input))
{

}

void FixedSourceSolver :: solve(){
}

}


