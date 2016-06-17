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

#include "solver_factory.hpp"

#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "eigen_solver.hpp"
#include "fixed_source_solver.hpp"
#include "monte_carlo_eigenvalue_solver.hpp"

namespace mocc {
SP_Solver_t SolverFactory(const pugi::xml_node &input, const CoreMesh &mesh)
{
    LogFile << "Initializing solver..." << std::endl;

    SP_Solver_t solver;

    if (input.empty()) {
        throw EXCEPT("No input specified for the solver.");
    }
    std::string type = input.attribute("type").value();
    if (type == "eigenvalue") {
        solver = std::make_shared<EigenSolver>(input, mesh);
    }
    else if (type == "fixed_source") {
        solver = std::make_shared<FixedSourceSolver>(input, mesh);
    }
    else if (type == "eigenvalue_mc") {
        solver = std::make_shared<mc::MonteCarloEigenvalueSolver>(input, mesh);
    }
    else {
        throw EXCEPT("Unrecognized solver type.");
    }

    LogFile << "Done initializing solver." << std::endl;

    return solver;
}
}
