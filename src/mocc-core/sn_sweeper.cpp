#include "sn_sweeper.hpp"

#include "error.hpp"

using std::cout;
using std::endl;

namespace mocc {
    SnSweeper::SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh,
            SP_XSMesh_t xs_mesh):
        xs_mesh_hom_( mesh ),
        mesh_hom_( mesh.n_pin(), 
                   mesh.n_pin(), 
                   mesh.nx(), 
                   mesh.ny(), 
                   mesh.nz() ),
        TransportSweeper( mesh_hom_, xs_mesh_ ),
        core_mesh_( &mesh ),
        ang_quad_( input.child("ang_quad") )
    {
        std::cout << "Generating Sn sweeper" << std::endl;
        
        // Make sure we have input from the XML
        if( input.empty() ) {
            EXCEPT("No input specified to initialize Sn sweeper.");
        }

        // Parse the number of inner iterations
        int int_in = input.attribute("n_inner").as_int(-1);
        if(int_in < 0) {
            EXCEPT("Invalid number of inner iterations specified (n_inner).");
        }
        n_inner_ = int_in;

        return;
    }

    void SnSweeper::sweep( int group ) {
        cout << __func__ << endl;
        return;
    }

    void SnSweeper::initialize() {
        return;
    }


    void SnSweeper::get_pin_flux( int ig, VecF& flux ) const { 
        return;
    }
}
