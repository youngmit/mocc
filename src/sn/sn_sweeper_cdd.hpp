#pragma once

#include <memory>

#include "pugixml.hpp"

#include "core_mesh.hpp"
#include "correction_data.hpp"
#include "sn_sweeper.hpp"

namespace mocc {
    class SnSweeper_CDD: public SnSweeper {
    public:
        SnSweeper_CDD( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        void set_corrections( CorrectionData *data ) {
            corrections_ = data;
        }
        
    private:
        template <typename CurrentWorker>
        void sweep_std( int group );

        const CorrectionData *corrections_;

        std::unique_ptr<const CorrectionData> my_corrections_;
    };
}
