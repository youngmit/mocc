#pragma once

#include "pugixml.hpp"

#include "global_config.hpp"
#include "transport_sweeper.hpp"
#include "angular_quadrature.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh, 
                SP_XSMesh_t xs_mesh );

        ~SnSweeper() { }

        void sweep( int group );

        void initialize();
        void get_pin_flux( int ig, VecF& flux ) const;

        void calc_fission_source( float_t k, 
                ArrayX& fission_source) const { }

        void output( H5File& file ) const { }
    private:
        XSMeshHomogenized xs_mesh_hom_;
        const CoreMesh* core_mesh_;
        Mesh mesh_hom_;
        unsigned int n_inner_;
        AngularQuadrature ang_quad_;
    };
}
