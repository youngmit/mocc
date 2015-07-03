#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "solver.hpp"
#include "core_mesh.hpp"
#include "transport_sweeper.hpp"

namespace mocc{

class FixedSourceSolver: public Solver{
public:
    // Initialize a FSS using the XML document
    FixedSourceSolver( const pugi::xml_node &input, const CoreMesh &mesh );

    ~FixedSourceSolver() {
        std::cout << "destroying fixed source solver" << std::endl;
    }
    
    // Solve a fixed source problem subject to the configuration in the XML
    // input. This can either be to some sort of tolerance, or for a single
    // group sweep
    void solve();
    
    // Set the group-independent fission source. The group-dependent fission
    // source is calculated using chi inside the group sweep.
    void set_FissionSource();

    // return the number of flat source regions
    int get_nReg(){
        return sweeper_->n_reg();
    }

private:
    UP_Sweeper_t sweeper_;
};
}
