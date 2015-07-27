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
        typedef std::vector< // group
                std::vector< // plane
                std::vector< // angle
                std::vector< float_t > > > > BCSet_t;   // BCs
        typedef std::vector< // plane
                std::vector< // angle
                std::vector< float_t > > > BCSet_Out_t; // BCs
    public:
        MoCSweeper( const pugi::xml_node &input,
                    const CoreMesh &mesh );
        
        ~MoCSweeper() { }
        
        void sweep(int group);

        void initialize();

        void calc_fission_source( float_t k,  ArrayX& fission_source) const;

    private:
        AngularQuadrature ang_quad_;
        RayData rays_;
        
        // Boundary condition. ordered by energy, plane, angle, ray
        BCSet_t boundary_;
        BCSet_Out_t boundary_out_;

        void sweep1g( int group );

        // Array of one group transport cross sections
        ArrayX xstr_;

        // Temporary storage for 1-group scalar flux
        ArrayX flux_1g_;

        // One-group, isotropic source, scaled by transport cross section
        ArrayX qbar_;

        // Number of inner iterations per group sweep
        unsigned int n_inner_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;




        // Update the boundary conditions
        void update_boundary( int group );
    };
}
