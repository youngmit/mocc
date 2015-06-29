#pragma once

#include <map>

#include "pugixml.hpp"

#include "pin_mesh.hpp"
#include "pin.hpp"
#include "material_lib.hpp"

namespace mocc {
    // The core mesh stores everything needed to represent the physical state
    // of the system. Pin meshes, material library, actual pin types, lattices
    // etc. The CoreMesh is then used to perform complex operations like ray
    // tracing, generation of coarse mesh, etc. A lot of the heavy lifting for
    // input processing happens in the constructor, and the CoreMesh assumes 
    // ownership of a lot of the structures used to represent the system.
    class CoreMesh {
    public:
        CoreMesh(pugi::xml_node &input);
    private:
        // Map for storing pin mesh objects indexed by user-specified IDs
        std::map<const int, PinMesh*> pin_meshes_;
        // The material library
        MaterialLib mat_lib_;
        // Map of actual pin objects, indexed by user-specified IDs
        std::map<const int, Pin*> pins_;
    };
}
