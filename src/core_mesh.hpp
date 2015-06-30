#pragma once

#include <map>
#include <memory>

#include "pugixml.hpp"

#include "pin_mesh.hpp"
#include "pin.hpp"
#include "material_lib.hpp"
#include "lattice.hpp"

namespace mocc {

    // The core mesh stores everything needed to represent the physical state
    // of the system. Pin meshes, material library, actual pin types, lattices
    // etc. The CoreMesh is then used to perform complex operations like ray
    // tracing, generation of coarse mesh, etc. A lot of the heavy lifting for
    // input processing happens in the constructor, and the CoreMesh assumes 
    // ownership of a lot of the structures used to represent the system.
    class CoreMesh {
    public:
        CoreMesh() {
            // do nothing
            return;
        }

        // Construct a CoreMesh from XML input. This routine is responsible for
        // parsing many of the tags in the XML document: <mesh>, <pin>,
        // <material_lib>, <lattice>, <core>
        CoreMesh( pugi::xml_node &input );

        ~CoreMesh();

    private:
        // Map for storing pin mesh objects indexed by user-specified IDs
        std::map<int, UP_PinMesh_t> pin_meshes_;

        // The material library
        MaterialLib mat_lib_;

        // Map of actual pin objects, indexed by user-specified IDs
        std::map<int, UP_Pin_t> pins_;

        // Map of lattice objects
        std::map<int, Lattice> lattices_;
    };

    typedef std::shared_ptr<CoreMesh> SP_CoreMesh_t;
    typedef std::unique_ptr<CoreMesh> UP_CoreMesh_t;
}
