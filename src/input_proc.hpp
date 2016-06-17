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

#include <map>
#include <memory>
#include <string>

#include "pugixml.hpp"

#include "core/core_mesh.hpp"
#include "core/timers.hpp"

#include "solvers/solver_factory.hpp"

namespace mocc{
    /**
    * InputProc is essentially a global storage location for the CoreMesh and
    * top-level Solver. Following construction, the "driver" can extract the
    * top-level Solver and call its Solver::solve() method.
    */
    class InputProcessor {
    public:
        /**
        * Given the filename of an XML document, parses the document into a tree
        * structure, then uses it to generate a CoreMesh and top-level Solver.
        */
        InputProcessor( std::vector<std::string> args );

        /**
         * \brief Actually process the contents of the file and construct
         * associated objects.
         */
        void process();

        /**
        * Return a shared pointer to the CoreMesh.
        */
        SP_CoreMesh_t core_mesh() {
            return core_mesh_;
        }

        /**
        * Return a shared pointer to the top-level solver.
        */
        SP_Solver_t solver() {
            return solver_;
        }

        /**
         * \brief Return the case name
         *
         * The case name assumes a value of the input file name, without the
         * '.xml' extension, unless another case name is specified using a
         * \<case_name\> tag in the input file.
         */
        std::string case_name() const {
            return case_name_;
        }

        const pugi::xml_document &document() const {
            return doc_;
        }

    private:
        // Timer for all input processing activities
        Timer &timer_;

        // Master core mesh object. Can be passed back to the driver or wherever
        SP_CoreMesh_t core_mesh_;

        // Top-level solver
        SP_Solver_t solver_;

        // XML document
        pugi::xml_document doc_;

        std::vector<std::string> args_;

        std::string case_name_;

    };
}
