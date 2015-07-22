#pragma once

#include <memory>

#include "pugixml.hpp"

#include "pin_mesh_base.hpp"
#include "pin_mesh_cyl.hpp"
#include "pin_mesh_rect.hpp"

namespace mocc {
    typedef std::shared_ptr<PinMesh> SP_PinMesh_t;
    typedef std::unique_ptr<PinMesh> UP_PinMesh_t;

    PinMesh* PinMeshFactory(const pugi::xml_node &input);
}
