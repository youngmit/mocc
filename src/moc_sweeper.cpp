#include "moc_sweeper.hpp"

#include "error.hpp"

namespace mocc {
    
    namespace {
        void sweep_moc(int group) {
            return;
        }
    }
    
    MoCSweeper :: MoCSweeper( const pugi::xml_node &input, 
                              const CoreMesh &mesh ):
        ang_quad_( input.child("ang_quad") ),
        rays_( input.child("rays"), ang_quad_, mesh )
    {
        
        // Make sure we have input from the XML
        if (input.empty()) {
            Error("No input specified to initialize MoC sweeper.");
        }

        return;
    }

    void MoCSweeper :: sweep(int group) {
        sweep_moc(group);
        return;
    }
}
