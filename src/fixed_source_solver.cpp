#include "fixed_source_solver.hpp"

#include "transport_sweeper_factory.hpp"

namespace mocc{
FixedSourceSolver :: FixedSourceSolver( const pugi::xml_node &input, 
        const CoreMesh &mesh )
{
    std::cout << "creating fixed source solver" << std::endl;
    std::cout << "transport sweeper located at " << sweeper_.get() << std::endl;
    sweeper_ = UP_Sweeper_t((std::move(TransportSweeperFactory(input, mesh))));
    std::cout << "transport sweeper located at " << sweeper_.get() << std::endl;
}

void FixedSourceSolver :: solve(){
}

}


