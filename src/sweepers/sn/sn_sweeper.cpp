#include "sn_sweeper.hpp"

#include <algorithm>
#include <string>

#include "pugixml.hpp"

#include "files.hpp"
#include "string_utils.hpp"

namespace mocc {
namespace sn {
    SnSweeper::SnSweeper( const pugi::xml_node &input, const CoreMesh& mesh ):
            TransportSweeper( input ),
            mesh_( mesh ),
            bc_type_( mesh.boundary() ),
            flux_1g_( ),
            xstr_( mesh.n_pin() ),
            q_( mesh_.n_pin() ),
            bc_in_( mesh.mat_lib().n_group(), ang_quad_, mesh_ ),
            bc_out_( 1, ang_quad_, mesh_ ),
            gs_boundary_( true )
        {
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
            }
            return;
        }

    void SnSweeper::output( H5::CommonFG *node ) const {
        auto dims = mesh_.dimensions();
        std::reverse( dims.begin(), dims.end() );

        // Make a group in the file to store the flux
        node->createGroup("flux");

        ArrayB2 flux = this->get_pin_flux();
        Normalize( flux.begin(), flux.end() );

        for( unsigned int ig=0; ig<n_group_; ig++ ) {
            std::stringstream setname;
            setname << "flux/" << std::setfill('0') << std::setw(3) << ig+1;

            ArrayB1 flux_1g = flux(blitz::Range::all(), ig);

            HDF::Write( node, setname.str(), flux_1g.begin(), flux_1g.end(),
                    dims);
        }

        LogFile << "Sn Sweeper:" << std::endl;

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
