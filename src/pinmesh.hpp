#pragma once

#include "pinmesh_base.hpp"
#include "pinmesh_cyl.hpp"
#include "pinmesh_rect.hpp"
#include "pugixml.hpp"
#include <memory>


typedef std::shared_ptr<PinMesh> SP_PinMesh;

SP_PinMesh PinMeshFactory(const pugi::xml_node &input);