#include "sn_sweeper.hpp"

#include <algorithm>
#include <string>

#include "pugixml.hpp"

#include "files.hpp"
#include "string_utils.hpp"


namespace {
    using namespace mocc;
    BC_Size_t boundary_helper( const Mesh &mesh ) {
        BC_Size_t bc_size = { (int)mesh.ny()*(int)mesh.nz(),
                              (int)mesh.nx()*(int)mesh.nz(),
                              (int)mesh.nx()*(int)mesh.ny() };
        return bc_size;
    }
}


namespace mocc {
namespace sn {
    SnSweeper::SnSweeper( const pugi::xml_node &input, const CoreMesh& mesh ):
            TransportSweeper( input ),
            timer_( RootTimer.new_timer("Sn Sweeper", true) ),
            timer_init_( timer_.new_timer("Initialization", true) ),
            timer_sweep_( timer_.new_timer("Sweep", false) ),
            mesh_( mesh ),
            bc_type_( mesh.boundary() ),
            flux_1g_( ),
            xstr_( mesh.n_pin() ),
            bc_in_( mesh.mat_lib().n_group(), ang_quad_, bc_type_,
                    boundary_helper(mesh) ),
            bc_out_( 1, ang_quad_, bc_type_, boundary_helper(mesh) ),
            gs_boundary_( true )
    {
        LogFile << "Constructing a base Sn sweeper" << std::endl;

        // Set up the cross-section mesh. If there is <data> specified, try to
        // use that, otherwise generate volume-weighted cross sections
        if( input.child("data").empty() ) {
            xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh) );
        } else {
            try {
                xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh, input) );
            }
            catch( Exception e ) {
                std::cerr << e.what() << std::endl;
                throw EXCEPT("Failed to create XSMesh for Sn Sweeper.");
            }
        }

        core_mesh_ = &mesh;
        n_reg_ = mesh.n_pin();
        n_group_ = xs_mesh_->n_group();
        flux_.resize( n_reg_, n_group_ );
        flux_old_.resize( n_reg_, n_group_ );
        vol_.resize( n_reg_ );

        // Set the mesh volumes. Same as the pin volumes
        int ipin = 0;
        for( auto &pin: mesh_ ) {
            int i = mesh_.coarse_cell( mesh_.pin_position(ipin) );
            vol_[i] = pin->vol();
            ipin++;
        }

        // Make sure we have input from the XML
        if( input.empty() ) {
            throw EXCEPT("No input specified to initialize Sn sweeper.");
        }

        // Parse the number of inner iterations
        int int_in = input.attribute("n_inner").as_int(-1);
        if( int_in < 0 ) {
            throw EXCEPT("Invalid number of inner iterations specified "
                    "(n_inner).");
        }
        n_inner_ = int_in;

        // Try to read boundary update option
        if( !input.attribute("boundary_update").empty() ) {
            std::string in_string =
                input.attribute("boundary_update").value();
            sanitize(in_string);

            if( (in_string == "gs") || (in_string == "gauss-seidel") ) {
                gs_boundary_ = true;
            } else if ( (in_string == "jacobi") || (in_string == "j") ) {
                gs_boundary_ = false;
            } else {
                throw EXCEPT("Unrecognized option for BC update!");
            }

            // For now, the BC doesnt support parallel boundary updates, so
            // disable Gauss-Seidel boundary update if we are using multiple
            // threads.
            /// \todo Add support for multi-threaded G-S boundary update in Sn
            if( (omp_get_max_threads() > 1) && gs_boundary_ ) {
                gs_boundary_ = false;
                LogScreen << "WARNING: Disabling Gauss-Seidel boundary update "
                    "in parallel Sn" << std::endl;
            }
        }

        timer_.toc();
        timer_init_.toc();

        return;
    }

    void SnSweeper::output( H5Node &node ) const {
        auto dims = mesh_.dimensions();
        std::reverse( dims.begin(), dims.end() );

        // Make a group in the file to store the flux
        node.create_group("flux");

        ArrayB2 flux = this->get_pin_flux();
        Normalize( flux.begin(), flux.end() );

        for( unsigned int ig=0; ig<n_group_; ig++ ) {
            std::stringstream setname;
            setname << "flux/" << std::setfill('0') << std::setw(3) << ig+1;

            ArrayB1 flux_1g = flux(blitz::Range::all(), ig);

            node.write( setname.str(), flux_1g.begin(), flux_1g.end(),
                    dims);
        }

        LogFile << "Sn Sweeper:" << std::endl;
        LogFile << "Angular Quadrature:" << std::endl;
        LogFile << ang_quad_ << std::endl;

        LogFile << "Boundary update: ";
        if( gs_boundary_ ) {
            LogFile << "Gauss-Seidel" << std::endl;
        } else {
            LogFile << "Jacobi" << std::endl;
        }

        LogFile << std::endl;

        xs_mesh_->output( node );
        return;
    }
} // namespace sn
} // namespace moc
