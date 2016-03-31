/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include "pugixml.hpp"

#include "solver.hpp"

namespace mocc {
    /**
     * In a more complicated world, this would interrogate the input XML to
     * determine the type of highest-level \ref Solver to use, allocate and
     * construct that \ref Solver and return a shared pointer to it. In
     * practice, we are only making instances of \ref EigenSolver.
     */
    SP_Solver_t SolverFactory( const pugi::xml_node &input,
        const CoreMesh &mesh );
}
