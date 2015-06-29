#pragma once

#include "solver.hpp"
#include "pugixml.hpp"
#include "fixedsrc_solver.hpp"

namespace mocc{

class EigenSolver: public Solver{
public:
    EigenSolver(pugi::xml_node input);
    void solve();

private:
    FixedSourceSolver m_fss;
};
};
