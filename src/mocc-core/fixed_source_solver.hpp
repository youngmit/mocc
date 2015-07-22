#pragma once

#include "pugixml.hpp"

#include "solver.hpp"
#include "source.hpp"
#include "core_mesh.hpp"
#include "transport_sweeper.hpp"

namespace mocc{
    class FixedSourceSolver: public Solver{
    public:
        // Initialize a FSS using the XML document
        FixedSourceSolver( const pugi::xml_node &input, const CoreMesh &mesh );
    
        ~FixedSourceSolver() {
        }
        
        // Solve a fixed source problem subject to the configuration in the XML
        // input. This can either be to some sort of tolerance, or for a single
        // group sweep
        void solve();

        void step();

        // Initialize the state of the FSS to start a new problem. For now this
        // just calls the same routine on the TransporSweeper, which in turn
        // initializes the scalar flux, boundary conditions, etc to some sort of
        // halfway-reasonable starting values.
        void initialize() {
            sweeper_->initialize();
        }
        
        // Set the group-independent fission source. The group-dependent fission
        // source is calculated internally
        void set_fission_source( const MatrixX* fs) {
std::cout << "Assigning fission source." << std::endl;
            fs_ = fs;
        }
    
        // return the number of flat source regions
        unsigned int n_reg() {
            return sweeper_->n_reg();
        }

        // Return a constant reference to the transport sweeper
        const TransportSweeper* sweeper() const {
            return sweeper_.get();
        }
    
    private:
        UP_Sweeper_t sweeper_;
        Source source_;
        // Pointer to the group-independent fission source. Usually comes from
        // an eigenvalue solver, if present
        const MatrixX* fs_;
        unsigned int ng_;
    };
}
