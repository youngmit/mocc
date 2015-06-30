#pragma once

#include "pin_mesh_base.hpp"
#include "pugixml.hpp"

namespace mocc {
    class PinMesh_Rect: public PinMesh{
    public:
    	PinMesh_Rect(const pugi::xml_node &input);
    private:
        VecF hx_;
        VecF hy_;
    };
}
