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
        
        // Set the group-independent fission source. The group-dependent fission
        // source is calculated internally
        void set_fission_source( const MatrixX* fs) {
            fs_ = fs;
        }
    
        // return the number of flat source regions
        unsigned int n_reg() {
            return sweeper_->n_reg();
        }
    
    private:
        unsigned int ng_;
        UP_Sweeper_t sweeper_;
        Source source_;
        // Pointer to the group-independent fission source. Usually comes from
        // an eigenvalue solver, if present
        const MatrixX* fs_;
    };
}
