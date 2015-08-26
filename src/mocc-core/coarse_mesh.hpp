#pragma once

#include <string>
#include <map>

#include "global_config.hpp"
#include "mesh.hpp"
#include "core_mesh.hpp"

namespace mocc {
    // The CoarseMesh type is sort of a one-stop shop for the data that might be
    // necessary for solvers that operate on pin-homogenized meshes. It doesnt
    // really store that much in the way of mesh properties since we guarantee
    // that it conforms to the CoreMesh, and is therefore a simple structured
    // grid. At some point it might be nice to extend the base Mesh type and use
    // this as the Sn mesh, too.
    class CoarseMesh {
    public:
        CoarseMesh(){ return; };
    private:
        // The actual data. Since what needs to be stored is problem-dependent,
        // we use a map, keyed on strings to keep things dynamic.
        std::map<std::string, VecF> data_;
    };
}
