#pragma once

#include "pugixml.hpp"

#include "global_config.hpp"
#include "moc_sweeper_2d3d.hpp"
#include "sn_sweeper_cdd.hpp"

namespace mocc {
    /**
     * This is an implementation of the 2D3D method. Each plane is treated with
     * a 2-D MoC sweeper, which produces the correction factors needed to treat
     * the entire system with a 3-D corrected diamond difference Sn sweeper.
     */
    class PlaneSweeper_2D3D: public TransportSweeper {
    public:
        PlaneSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        void initialize();

        void get_pin_flux( int ig, VecF &flux ) const;

        void calc_fission_source( float_t k, ArrayX &fission_source ) const;
    
        void output( H5File& file ) const;

        void homogenize( CoarseData &data ) const {
            
        }

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return sn_sweeper_.get_homogenized_xsmesh();
        }
        
    private:
        SnSweeper_CDD sn_sweeper_;
        MoCSweeper_2D3D moc_sweeper_;
    };
}
