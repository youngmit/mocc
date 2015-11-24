#include "sn_sweeper.hpp"

#include <iomanip>

#include "mocc-core/error.hpp"
#include "mocc-core/files.hpp"
#include "mocc-core/utils.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
    SnSweeper::SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh ):
        mesh_( mesh ),
        ang_quad_( input.child("ang_quad") ),
        bc_type_( mesh.boundary() ),
        flux_1g_( mesh_.n_pin() ),
        xstr_( mesh.n_pin() ),
        q_( mesh_.n_pin() ),
        bc_in_( mesh.mat_lib().n_group(), ang_quad_, mesh_ ),
        bc_out_( 1, ang_quad_, mesh_ )
    {
        LogFile << "Constructing a base Sn sweeper" << std::endl;

        // Set up all of the stuff that would normally be done by the
        // TransportSweeper constructor. There is probably a better and more
        // maintainable way to do this; will revisit.
        core_mesh_ = &mesh;
        xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh) );
        n_reg_ = mesh.n_pin();
        n_group_ = xs_mesh_->n_group();
        flux_.resize( n_reg_, n_group_ );
        flux_old_.resize( n_reg_, n_group_ );
        vol_.resize( n_reg_ );

        // Set the mesh volumes. Same as the pin volumes
        int ipin = 0;
        for( auto &pin: mesh_ ) {
            int i = mesh_.index_lex( mesh_.pin_position(ipin) );
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

        return;
    }

    void SnSweeper::homogenize( CoarseData &data ) const {
        return;
    }

    void SnSweeper::initialize() {
        flux_ = 1.0;
        flux_old_ = 1.0;
        bc_in_.initialize(1.0/FPI);

        return;
    }

    void SnSweeper::get_pin_flux_1g( int ig, VecF& flux ) const {
        size_t n = mesh_.n_pin();
        flux.clear();

        auto tmp_f = flux_(blitz::Range::all(), ig);
        for( const auto &f: tmp_f ) {
            flux.push_back(f);
        }

        return;
    }

    void SnSweeper::output( H5::CommonFG *node ) const {
        auto dims = mesh_.dimensions();
        std::reverse( dims.begin(), dims.end() );

        // Make a group in the file to store the flux
        node->createGroup("flux");

        VecF flux = this->get_pin_flux();
        Normalize( flux.begin(), flux.end() );

        auto flux_it = flux.cbegin();

        for( unsigned int ig=0; ig<n_group_; ig++ ) {
            std::stringstream setname;
            setname << "flux/" << std::setfill('0') << std::setw(3) << ig+1;

            flux_it = HDF::Write( node, setname.str(), flux_it,
                    flux_it+mesh_.n_pin(), dims);
        }

        xs_mesh_->output( node );
        return;
    }

    void SnSweeper::check_balance( int group ) const {
        if( !coarse_data_ ) {
            throw EXCEPT("No coarse data. Need it to look at currents.");
        }
        for( size_t icell=0; icell<mesh_.n_pin(); icell++ ) {
            real_t b = 0.0;

            // Current
            b -= coarse_data_->current(
                    mesh_.coarse_surf(icell, Surface::EAST), group ) * 
                    mesh_.coarse_area( icell, Surface::EAST );
            b -= coarse_data_->current(
                    mesh_.coarse_surf(icell, Surface::NORTH), group ) *
                    mesh_.coarse_area( icell, Surface::NORTH );
            b -= coarse_data_->current(
                    mesh_.coarse_surf(icell, Surface::TOP), group ) *
                    mesh_.coarse_area( icell, Surface::TOP );
            b += coarse_data_->current(
                    mesh_.coarse_surf(icell, Surface::WEST), group ) *
                    mesh_.coarse_area( icell, Surface::WEST );
            b += coarse_data_->current(
                    mesh_.coarse_surf(icell, Surface::SOUTH), group ) *
                    mesh_.coarse_area( icell, Surface::SOUTH );
            b += coarse_data_->current(
                    mesh_.coarse_surf(icell, Surface::BOTTOM), group ) *
                    mesh_.coarse_area( icell, Surface::BOTTOM );            

            // Source
            b += (*source_)[icell]*vol_[icell];

            // Internal removal
            b -= flux_1g_(icell) * 
                (*xs_mesh_)[icell].xsmacrm()[group] * vol_[icell];

            cout << "Cell balance: " << b << endl;
            cout << endl;
        }
    }
}
