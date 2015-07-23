#include "fixed_source_solver.hpp"

#include "error.hpp"
#include "transport_sweeper_factory.hpp"

namespace mocc {
    FixedSourceSolver::FixedSourceSolver( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        sweeper_( UP_Sweeper_t( TransportSweeperFactory(input, mesh) ) ),
        source_( mesh.n_reg(), sweeper_->xs_mesh(), sweeper_->cflux() ),
        fs_( nullptr ),
        ng_( sweeper_->n_grp() )
    {
        sweeper_->assign_source( &source_ );
        return;
    }

    // Perform source iteration
    void FixedSourceSolver::solve() {
        // not implemented yet
        Error("Stand-alone source iteration is not implemented yet.");
    }

    // Perform a single group sweep
    void FixedSourceSolver::step() {
        for( unsigned int ig=0; ig<ng_; ig++ ) {
std::cout << "sweeping group: " << ig << std::endl;
            // Set up the source
            if( fs_ == nullptr ) {
                Error("No fission source associated!");
            }
            source_.fission( *fs_, ig );
            source_.in_scatter( ig );

            sweeper_->sweep(ig);
        }
    }
}
