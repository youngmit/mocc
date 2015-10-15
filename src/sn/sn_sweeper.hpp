#pragma once

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "coarse_data.hpp"
#include "core_mesh.hpp"
#include "global_config.hpp"
#include "mesh.hpp"
#include "sn_boundary.hpp"
#include "sn_current_worker.hpp"
#include "sn_source.hpp"
#include "transport_sweeper.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh );

        ~SnSweeper() { }

        virtual void sweep( int group );

        void initialize();

        void get_pin_flux_1g( int ig, VecF& flux ) const;

        void output( H5::CommonFG *node ) const;

        // Override the create_source() method to make an SnSource instead of
        // the regular
        UP_Source_t create_source() const {
            Source *s = new SnSource( n_reg_, xs_mesh_.get(), this->cflux());
            UP_Source_t source( s );
            return source;
        }

        void homogenize( CoarseData &data ) const;

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return std::static_pointer_cast<XSMeshHomogenized>( xs_mesh_ );
        }

    protected:
        const CoreMesh &mesh_;

        // Mesh parameters
        int nx_;
        int ny_;
        int nz_;
        VecF hx_;
        VecF hy_;
        VecF hz_;

        // Update the boundary conditions 
        void update_boundary( int group );

        unsigned int n_inner_;
        AngularQuadrature ang_quad_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;
        
        // Temporary storage for 1-group scalar flux
        ArrayX flux_1g_;
        
        // Temporary storage of the current-group transport cross section
        ArrayX xstr_;

        // Single-group isotropic source, should include in-scatter
        ArrayX q_;

        // Incomming boundary condition
        SnBoundary bc_in_;

        // Outgoing boundary condition. Only difined for one group
        SnBoundary bc_out_;

        template <typename CurrentWorker>
        void sweep_dd( int group );
    };

    typedef std::unique_ptr<SnSweeper> UP_SnSweeper_t;
}
