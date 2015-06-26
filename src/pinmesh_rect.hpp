#pragma once

#include "pinmesh_base.hpp"
#include "pugixml.hpp"

class PinMesh_Rect: public PinMesh{
public:
	PinMesh_Rect(const pugi::xml_node &input);
private:
};