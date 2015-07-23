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
            float_t fw_start;
            float_t bw_start;
            float_t fw_end;
            float_t bw_end;
        };
        typedef std::vector< 
                std::vector<
                std::vector<
                std::vector< BC > > > > BCSet_t;
        typedef std::vector< std::vector< BC > > BCRays_t;
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

        // Array of one group transport cross sections
        ArrayX xstr_;

        // Temporary storage for 1-group scalar flux
        ArrayX flux_1g_;

        // One-group, isotropic source, scaled by transport cross section
        ArrayX qbar_;

        // Number of inner iterations per group sweep
        unsigned int n_inner_;
    };
}
