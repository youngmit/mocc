#pragma once

#include "correction_data.hpp"
#include "global_config.hpp"
#include "moc_sweeper.hpp"

namespace mocc {
    class MoCSweeper_2D3D: public MoCSweeper {
    public:
        MoCSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        void set_coupling( CorrectionData *data, 
                const XSMeshHomogenized *xsmesh) {
            corrections_ = data;
            sn_xs_mesh_ = xsmesh;
        }

    private:
        void sweep1g_final( int group );

        /**
         * Given homogenized angular flux and total cross section data,
         * calculate the correction factors for CDD
         */
        void calculate_corrections( size_t ang, size_t group, ArrayF flux_surf,
                ArrayF flux_node, ArrayF sigt );

        CorrectionData* corrections_;

        const XSMeshHomogenized* sn_xs_mesh_;
    };
}
