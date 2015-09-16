#pragma once

#include "pugixml.hpp"

#include "core_mesh.hpp"
#include "correction_data.hpp"
#include "sn_sweeper.hpp"

namespace mocc {
    class SnSweeper_CDD: public SnSweeper {
    public:
        SnSweeper_CDD( const pugi::xml_node &input, const CoreMesh &mesh );

        void set_corrections( CorrectionData *data ) {
            corrections_ = data;
        }
        
    private:
        const CorrectionData *corrections_;
    };
}
