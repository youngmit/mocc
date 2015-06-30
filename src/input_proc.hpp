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
        SP_CoreMesh_t core_mesh_;
    };
}
