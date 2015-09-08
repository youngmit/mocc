#pragma once

#include "moc_sweeper.hpp"
namespace mocc {
    class MoCSweeper_2D3D: public MoCSweeper {
    public:
        MoCSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );
    };
}
