#pragma once

#include "mocc-core/global_config.hpp"

#include "sweepers/moc/moc_sweeper.hpp"

#include "sweepers/sn/correction_data.hpp"

namespace mocc {
    class MoCSweeper_2D3D: public moc::MoCSweeper {
    public:
        MoCSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        void set_coupling( std::shared_ptr<CorrectionData> data,
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

        std::shared_ptr<CorrectionData> corrections_;

        const XSMeshHomogenized* sn_xs_mesh_;
    };
}
