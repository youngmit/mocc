#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "eigen_interface.hpp"
#include "transport_sweeper.hpp"
#include "ray_data.hpp"
#include "core_mesh.hpp"
#include "angular_quadrature.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    class MoCSweeper: public TransportSweeper{
        struct BC {
            float_t fw;
            float_t bw;
        };
        typedef std::vector< 
                std::vector<
                std::vector<
                std::vector< BC > > > > BCSet_t;
    public:
        MoCSweeper( const pugi::xml_node &input,
                    const CoreMesh &mesh );
        
        ~MoCSweeper() {
        }
        
        void sweep(int group);

        void initialize();

        void calc_fission_source( float_t k,  ArrayX& fission_source) const;
        
    private:
        AngularQuadrature ang_quad_;
        RayData rays_;
        
        // Boundary condition. ordered by energy, plane, angle, ray
        BCSet_t boundary_;

        void sweep1g( int group );
    };
}
