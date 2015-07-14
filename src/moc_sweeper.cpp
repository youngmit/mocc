#include "moc_sweeper.hpp"

#include "error.hpp"

namespace mocc {
    
    MoCSweeper :: MoCSweeper( const pugi::xml_node &input, 
                              const CoreMesh &mesh ):
        ang_quad_( input.child("ang_quad") ),
        rays_( input.child("rays"), ang_quad_, mesh )
    {   
        // Make sure we have input from the XML
        if( input.empty() ) {
            Error("No input specified to initialize MoC sweeper.");
        }

        xs_mesh_ = XSMesh(mesh);   

        return;
    }

    void MoCSweeper :: sweep(int group) {
        return;
    }

    void MoCSweeper::sweep1g( int group ) {


        return;
    }
}
