#include "fixed_source_solver.hpp"

#include "transport_sweeper_factory.hpp"

#include "core/error.hpp"
#include "core/files.hpp"
#include "core/h5file.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
    FixedSourceSolver::FixedSourceSolver( const pugi::xml_node &input,
            const CoreMesh &mesh ) try :
        sweeper_( UP_Sweeper_t( TransportSweeperFactory(input, mesh) ) ),
        source_( sweeper_->create_source(input.child("source")) ),
        fs_( nullptr ),
        ng_( sweeper_->n_group() ),
        fixed_source_( false )
    {
        std::string type = input.attribute("type").value();
        // See if we are creating a fully-specified FSS. If the passed-in input
        // is type="fixed_source" do extra stuff.
        if( type == "fixed_source" ) {
            fixed_source_ = true;
            LogFile << "Using an explicitly-defined fixed source solver"
                << std::endl;

            // Iteration limit
            max_iter_ = input.attribute("max_iter").as_int(-1);
            if( max_iter_ <= 0 ) {
                throw EXCEPT( "Failed to parse reasonable number of maximum "
                        "iterations." );
            }

            // Convergence criterion
            flux_tol_ = input.attribute("flux_tol").as_float(-1.0);
            if( flux_tol_ <= 0.0 ) {
                throw EXCEPT( "Failed to parse a reasonable flux tolerance." );
            }

            LogFile << "Maximum number of outer iterations: " << max_iter_ <<
                std::endl;
            LogFile << "Flux tolerance: " << flux_tol_ << std::endl;

            // Source
            if( input.child("source").empty() ) {
                throw EXCEPT("Top-level fixed source solver needs an explicit "
                        "source!");
            }
            source_->add_external( input.child("source") );
        }

        sweeper_->assign_source( source_.get() );
        return;
    }
    catch (Exception e) {
        Fail(e);
    }

    // Perform source iteration
    void FixedSourceSolver::solve() {
        for( size_t iouter=0; iouter<max_iter_; iouter++ ) {
            this->step();

            real_t resid = sweeper_->flux_residual();
            cout << iouter << " " << resid << endl;

            if( resid < flux_tol_ ) {
                break;
            }
        }
    }

    // Perform a single group sweep
    void FixedSourceSolver::step() {
        // Tell the sweeper to stash its old flux
        sweeper_->store_old_flux();
        for( size_t ig=0; ig<ng_; ig++ ) {
            // Set up the source
            source_->initialize_group( ig );
            if( fs_ ) {
                source_->fission( *fs_, ig );
            }

            source_->in_scatter( ig );

            sweeper_->sweep( ig );
        }
    }

    void FixedSourceSolver::output( H5Node &node ) const {
        // Provide energy group upper bounds
        // We do this here, to prevent collisions between possibly-multiple
        // sweepers colliding.
        node.write( "ng", (int)sweeper_->n_group() );
        node.write( "eubounds", sweeper_->xs_mesh().eubounds(), VecI(1, ng_) );
        sweeper_->output( node );
        return;
    }
}
