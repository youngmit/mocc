#pragma once

#include "pugixml.hpp"

#include "solver.hpp"
#include "source.hpp"
#include "core_mesh.hpp"
#include "transport_sweeper.hpp"
#include "h5file.hpp"

namespace mocc{
    /**
    * This Solver (FSS) attempts to solve the fixed source problem. For now, the
    * fixed source must be provided by some solver above the FSS, in the form of
    * a Source object, however in the future it might be useful to be able to
    * supply a user-defined Source for non-eigenvalue problems.
    * 
    * Right now, the FSS is used by the EigenSolver to converge the flux
    * solution for intermediate "fixed" sources for each eigenvalue step.
    */
    class FixedSourceSolver: public Solver {
    public:
        /**
        * Initialize a FSS using an XML node and CoreMesh. The expects the
        * passed XML node to be a valid \<solver\> tag containing a relevant
        * \<sweeper\> tag, which is needed by the TransportSweeperFactory() to
        * generate a TransportSweeper.
        */
        FixedSourceSolver( const pugi::xml_node &input, const CoreMesh &mesh );
    
        ~FixedSourceSolver() {
        }
        
        /**
        * For now, there is no actual implementation of this method, since there
        * is no functionality for specifying a user-defined Source. In practice,
        * the FSS is driven via the step() routine by the EigenSolver.
        *
        * Ideally, this would solve a fixed source problem subject to the
        * configuration in the XML input. This can either be to some sort of
        * tolerance, or for a fixed number of group sweeps.
        */
        void solve();

        /**
        * Instructs the sweeper to store the old value of the flux, then
        * performs a group sweep.
        */
        void step();

        /**
        * Initialize the state of the FSS to start a new problem. For now this
        * just calls the same routine on the TransportSweeper, which in turn
        * initializes the scalar flux, boundary conditions, etc. to some sort of
        * halfway-reasonable starting values.
        */
        void initialize() {
            sweeper_->initialize();
        }
        
        /// Set the group-independent fission source. The group-dependent fission
        /// source is calculated internally by the Source object.
        void set_fission_source( const ArrayX* fs) {
            fs_ = fs;
        }
    
        /// Return the number of flat source regions.
        unsigned int n_reg() {
            return sweeper_->n_reg();
        }

        /// Return the number of energy groups
        unsigned int n_group() {
            return ng_;
        }

        /// Return a constant reference to the transport sweeper.
        const TransportSweeper* sweeper() const {
            return sweeper_.get();
        }

        void output( H5File& file ) const {
            sweeper_->output( file );
            return;
        }
    
    private:
        UP_Sweeper_t sweeper_;
        UP_Source_t source_;
        // Pointer to the group-independent fission source. Usually comes from
        // an eigenvalue solver, if present
        const ArrayX* fs_;
        unsigned int ng_;
    };
}
