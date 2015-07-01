#pragma once

#include <map>
#include <memory>

#include "core_mesh.hpp"
#include "solver_factory.hpp"

namespace mocc{
    class InputProc {
    public:
    	InputProc(const char* filename);
        SP_CoreMesh_t core_mesh() {
            return core_mesh_;
        }

    private:
        // Master core mesh object. Can be passed back to the driver or wherever
        SP_CoreMesh_t core_mesh_;

        // Top-level solver
        SP_Solver_t solver_;

    };
}
