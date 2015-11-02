#include "sn_sweeper.hpp"

#include <iomanip>

#include "error.hpp"
#include "files.hpp"
#include "utils.hpp"

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
        q_( mesh.n_pin() )
    {
        LogFile << "Constructing a base Sn sweeper" << std::endl;

        // Set up all of the stuff that would normally be done by the
        // TransportSweeper constructor. There is probably a better and more
        // maintainable way to do this; will revisit.
        core_mesh_ = &mesh;
        xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh) );
        n_reg_ = mesh.n_pin();
        n_group_ = xs_mesh_->n_group();
        flux_.resize( n_reg_ * n_group_ );
        flux_old_.resize( n_reg_ * n_group_ );
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
            throw EXCEPT("Invalid number of inner iterations specified (n_inner).");
        }
        n_inner_ = int_in;

        bc_in_ = SnBoundary( n_group_, ang_quad_.ndir(), mesh_.nx(), 
                mesh_.ny(), mesh_.nz() );
        bc_out_ = SnBoundary( 1, ang_quad_.ndir(), mesh_.nx(), mesh_.ny(), 
                mesh_.nz() );

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

        ArrayF tmp_f = flux_[std::slice(n*ig, n, 1)];
        for( const auto &f: tmp_f ) {
            flux.push_back(f);
        }

        return;
    }

    void SnSweeper::update_boundary( int group ) {
        int iang = 0;

        for( auto &ang: ang_quad_ ) {
            Surface surf;
            Normal norm;
            int iang_refl;

            // X-normal surface
            norm = Normal::X_NORM;
            surf = (ang.ox > 0) ? Surface::WEST : Surface::EAST;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_type_[(int)surf] == Boundary::REFLECT ) {
                auto new_bc = bc_out_.get_face( 0, iang_refl, norm );
                bc_in_.set_face( group, iang, norm, new_bc );
            } else {
                bc_in_.zero_face(group, iang, norm);
            }
            // Y-normal surface
            norm = Normal::Y_NORM;
            surf = (ang.oy > 0) ? Surface::SOUTH : Surface::NORTH;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_type_[(int)surf] == Boundary::REFLECT ) {
                auto new_bc = bc_out_.get_face( 0, iang_refl, norm );
                bc_in_.set_face(group, iang, norm, new_bc);
            } else {
                bc_in_.zero_face(group, iang, norm);
            }
            // Z-normal surface
            norm = Normal::Z_NORM;
            surf = (ang.oz > 0) ? Surface::BOTTOM : Surface::TOP;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_type_[(int)surf] == Boundary::REFLECT ) {
                auto new_bc = bc_out_.get_face( 0, iang_refl, norm );
                bc_in_.set_face(group, iang, norm, new_bc );
            } else {
                bc_in_.zero_face(group, iang, norm);
            }

            iang++;
        }
        return;
    }

    void SnSweeper::update_boundary( int group, int iang ) {
        Angle ang = ang_quad_[iang];
        Surface surf;
        Normal norm;
        int iang_refl = 0;

        // X-normal surface
        norm = Normal::X_NORM;
        surf = (ang.ox > 0) ? Surface::EAST : Surface::WEST;
        iang_refl = ang_quad_.reflect( iang, norm );
        if( bc_type_[(int)surf] == Boundary::REFLECT ) {
            auto new_bc = bc_out_.get_face( 0, iang, norm );
            bc_in_.set_face( group, iang_refl, norm, new_bc );
        } else {
            bc_in_.zero_face(group, iang_refl, norm);
        }
        // Y-normal surface
        norm = Normal::Y_NORM;
        surf = (ang.oy > 0) ? Surface::NORTH : Surface::SOUTH;
        iang_refl = ang_quad_.reflect( iang, norm );
        if( bc_type_[(int)surf] == Boundary::REFLECT ) {
            auto new_bc = bc_out_.get_face( 0, iang, norm );
            bc_in_.set_face(group, iang_refl, norm, new_bc);
        } else {
            bc_in_.zero_face(group, iang_refl, norm);
        }
        // Z-normal surface
        norm = Normal::Z_NORM;
        surf = (ang.oz > 0) ? Surface::TOP : Surface::BOTTOM;
        iang_refl = ang_quad_.reflect( iang, norm );
        if( bc_type_[(int)surf] == Boundary::REFLECT ) {
            auto new_bc = bc_out_.get_face( 0, iang, norm );
            bc_in_.set_face(group, iang_refl, norm, new_bc );
        } else {
            bc_in_.zero_face(group, iang_refl, norm);
        }
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
}
