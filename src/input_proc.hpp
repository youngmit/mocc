#pragma once

#include <map>
#include <memory>

#include "mocc-core/core_mesh.hpp"
#include "mocc-core/solver_factory.hpp"

namespace mocc{
    /**
    * InputProc is essentially a global storage location for the CoreMesh and
    * top-level Solver. Following construction, the "driver" can extract the
    * top-level Solver and call its Solver::solve() method.
    */
    class InputProc {
    public:
        /**
        * Given the filename of an XML document, parses the document into a tree
        * structure, then uses it to generate a CoreMesh and top-level Solver.
        */
        InputProc(const char* filename);

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

    private:
        // Master core mesh object. Can be passed back to the driver or wherever
        SP_CoreMesh_t core_mesh_;

        // Top-level solver
        SP_Solver_t solver_;

    };
}
