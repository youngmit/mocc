#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "coarse_data.hpp"
#include "core_mesh.hpp"
#include "eigen_interface.hpp"
#include "ray_data.hpp"
#include "transport_sweeper.hpp"
#include "xs_mesh.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class MoCSweeper: public TransportSweeper {
        struct BC {
            real_t fw;
            real_t bw;
        };
        typedef std::vector< // group
                std::vector< // plane
                std::vector< // angle
                std::vector< real_t > > > > BCSet_t;   // BCs
        typedef std::vector< // plane
                std::vector< // angle
                std::vector< real_t > > > BCSet_Out_t; // BCs
    public:
        MoCSweeper( const pugi::xml_node &input,
                    const CoreMesh& mesh );
        
        ~MoCSweeper() { }
        
        virtual void sweep(int group);

        void initialize();

        void get_pin_flux( int group, VecF& flux ) const;

        void output( H5File& file ) const;

        void homogenize( CoarseData &data ) const;

        /**
         * Return a copy of the sweeper's angular quadrature.
         */
        AngularQuadrature get_ang_quad() const {
            return ang_quad_;
        }

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return SP_XSMeshHomogenized_t( 
                    new XSMeshHomogenized( mesh_ ) );
        }
    protected:
        const CoreMesh& mesh_;

        AngularQuadrature ang_quad_;
        RayData rays_;
        
        // Boundary condition. ordered by energy, plane, angle, ray
        BCSet_t boundary_;
        BCSet_Out_t boundary_out_;

        /// Perform a single "inner" iteration sweep
        void sweep1g( int group );

        /**
         * Perform a single "inner" iteration sweep, collecting special
         * information needed for coupling. At a bare minimum, this should
         * compute the interpin surface currents for use with CMFD, though
         * derived types of MoCSweeper may collect other data as needed.
         */
        virtual void sweep1g_final( int group );

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
