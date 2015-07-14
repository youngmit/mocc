#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "transport_sweeper.hpp"
#include "ray_data.hpp"
#include "core_mesh.hpp"
#include "angular_quadrature.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    class MoCSweeper: public TransportSweeper{
    public:
        MoCSweeper( const pugi::xml_node &input,
                    const CoreMesh &mesh );
        
        ~MoCSweeper() {
        }
        
        void sweep(int group);
        
    private:
        AngularQuadrature ang_quad_;
        RayData rays_;

        void sweep1g( int group );
    };
}
