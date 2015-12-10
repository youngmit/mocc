#pragma once

#include "pugixml.hpp"

#include "mocc-core/ray_data.hpp"
#include "mocc-core/core_mesh.hpp"

namespace mocc { namespace aux {
    void output_geometry( const pugi::xml_node &input, const CoreMesh &mesh );
} } //namespaces
