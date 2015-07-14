#include "fixed_source_solver.hpp"

#include "error.hpp"
#include "transport_sweeper_factory.hpp"

namespace mocc {
    FixedSourceSolver::FixedSourceSolver( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        sweeper_( UP_Sweeper_t(TransportSweeperFactory(input, mesh)) ),
        source_( mesh.n_reg(), sweeper_->xs_mesh() )
    {
        
    }

    // Perform source iteration
    void FixedSourceSolver::solve() {
        // not implemented yet
        Error("Stand-alone source iteration is not implemented yet.");
    }

    // Perform a single group sweep
    void FixedSourceSolver::step() {
        for( int ig=0; ig<ng_; ig++ ) {
            sweeper_->sweep(ig);
        }
    }
}


