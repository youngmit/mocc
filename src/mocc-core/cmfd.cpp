#include "cmfd.hpp"

namespace mocc {
    CMFD::CMFD( Mesh *mesh, int ng ):
        mesh_(mesh),
        coarse_data_( mesh->n_pin(), ng )
    {
        return;
    }
    
    
    void CMFD::solve() {
        

        return;
    }
}
