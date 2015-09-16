#include "cmfd.hpp"

namespace mocc {
    CMFD::CMFD( Mesh *mesh, SP_XSMeshHomogenized_t xsmesh):
        mesh_(mesh),
        xsmesh_( xsmesh ),
        coarse_data_( mesh->n_pin(), mesh->n_surf(), xsmesh->n_group() )
    {
        return;
    }
    
    
    void CMFD::solve() {
        
        // Construct the system matrix
        

        return;
    }
}
