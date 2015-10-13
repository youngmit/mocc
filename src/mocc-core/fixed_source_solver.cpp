#include "fixed_source_solver.hpp"

#include "error.hpp"
#include "transport_sweeper_factory.hpp"

namespace mocc {
    FixedSourceSolver::FixedSourceSolver( const pugi::xml_node &input, 
            const CoreMesh &mesh ) try :
        sweeper_( UP_Sweeper_t( TransportSweeperFactory(input, mesh) ) ),
        source_( sweeper_->create_source() ),
        fs_( nullptr ),
        ng_( sweeper_->n_group() )
    {
        
        sweeper_->assign_source( source_.get() );
        return;
    }
    catch (Exception e) {
            Fail(e);
    }

    // Perform source iteration
    void FixedSourceSolver::solve() {
        // not implemented yet
        Error("Stand-alone source iteration is not implemented yet.");
    }

    // Perform a single group sweep
    void FixedSourceSolver::step() {
        // Tell the sweeper to stash its old flux
        sweeper_->store_old_flux();
        for( unsigned int ig=0; ig<ng_; ig++ ) {
            // Set up the source
            if( fs_ == nullptr ) {
                Error("No fission source associated!");
            }
            source_->fission( *fs_, ig );
            source_->in_scatter( ig );

            sweeper_->sweep(ig);
        }
    }

    void FixedSourceSolver::output( H5::CommonFG *node ) const {
        // Provide energy group upper bounds
        // We do this here, to prevent collisions between possibly-multiple
        // sweepers colliding.
        HDF::Write( node, "ng", sweeper_->n_group() );
        HDF::Write( node, "eubounds", sweeper_->xs_mesh().eubounds(), 
                VecI(1, ng_) );
        sweeper_->output( node );
        return;
    }
}
