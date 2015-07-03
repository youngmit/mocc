#pragma once

#include "pugixml.hpp"

#include "solver.hpp"
#include "core_mesh.hpp"
#include "fixed_source_solver.hpp"

namespace mocc{

    class EigenSolver: public Solver{
    public:
        EigenSolver( const pugi::xml_node &input, const CoreMesh &mesh );
        void solve();
    
    private:
        FixedSourceSolver m_fss;
    };
}
