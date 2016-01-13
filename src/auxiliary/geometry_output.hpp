#pragma once

#include "pugixml.hpp"

#include "moc/ray_data.hpp"
#include "core/core_mesh.hpp"

namespace mocc { namespace aux {
    void output_geometry( const pugi::xml_node &input, const CoreMesh &mesh );
} } //namespaces
