#pragma once

#include "correction_data.hpp"
#include "moc_sweeper.hpp"

namespace mocc {
    class MoCSweeper_2D3D: public MoCSweeper {
    public:
        MoCSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );

                
        void set_corrections( CorrectionData *data ) {
            corrections_ = data;
        }

    private:
        /**
         * Override the sweep1g_final() routine on tha standard MoCSweeper to
         * also compute the corrected diamond difference correction factors.
         */
        void sweep1g_final( int group );

        CorrectionData* corrections_;
    };
}
