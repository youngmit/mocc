#pragma once

#include "pin_mesh_base.hpp"
#include "pin_mesh_cyl.hpp"
#include "pin_mesh_rect.hpp"
#include "pugixml.hpp"
#include <memory>


typedef std::shared_ptr<PinMesh> SP_PinMesh;

SP_PinMesh PinMeshFactory(const pugi::xml_node &input);