#include "sn_sweeper.hpp"

#include <iomanip>

#include "error.hpp"
#include "utils.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
    SnSweeper::SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh ):
        mesh_( mesh ),
        ang_quad_( input.child("ang_quad") ),
        bc_type_( mesh.boundary() ),
        xstr_( mesh.n_pin(), 1 ), 
        q_( mesh.n_pin(), 1 )
    {
        // Set up all of the stuff that would normally be done by the
        // TransportSweeper constructor. There is probably a better and more
        // maintainable way to do this; will revisit.
        core_mesh_ = &mesh;
        xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh) );
        n_reg_ = mesh.n_pin();
        n_group_ = xs_mesh_->n_group();
        flux_.resize( n_reg_, n_group_ );
        flux_old_.resize( n_reg_, n_group_ );
        vol_.resize( n_reg_, 1 );

        flux_1g_.resize( n_reg_, 1 );

        // Set the mesh volumes. Same as the pin volumes
        int ipin = 0;
        for( auto &pin: mesh_ ) {
            int i = mesh_.index_lex( mesh_.pin_position(ipin) );
            vol_(i) = pin->vol();
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

        // Set up the mesh dimensions
        hx_ = mesh_.pin_dx();
        hy_ = mesh_.pin_dy();
        hz_ = mesh_.hz();

        nx_ = hx_.size();
        ny_ = hy_.size();
        nz_ = hz_.size();

        bc_in_ = SnBoundary( n_group_, ang_quad_.ndir(), nx_, ny_, nz_ );
        bc_out_ = SnBoundary( 1, ang_quad_.ndir(), nx_, ny_, nz_ );

        return;
    }

    void SnSweeper::homogenize( CoarseData &data ) const {
        return;
    }

    void SnSweeper::sweep( int group ) {
        // Store the transport cross section somewhere useful
        for( auto &xsr: *xs_mesh_ ) {
            real_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_(ireg) = xstr;
            }
        }
        
        flux_1g_ = flux_.col( group );

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            // Set the source (add upscatter and divide by 4PI)
            source_->self_scatter( group, flux_1g_, q_ );

            if( inner == n_inner_-1 && coarse_data_ ) {
                // Wipe out the existing currents
                coarse_data_->current.col( group ) = 0.0;
                this->sweep_dd<sn::Current>( group );
            } else {
                this->sweep_dd<sn::NoCurrent>( group );
            }
        }
        flux_.col( group ) = flux_1g_;

        return;
    }

    template <typename CurrentWorker>
    void SnSweeper::sweep_dd( int group ) {
        CurrentWorker cw( coarse_data_, &mesh_ );
        flux_1g_.fill(0.0);

		ArrayF x_flux(ny_*nz_);
		ArrayF y_flux(nx_*nz_);
		ArrayF z_flux(nx_*ny_);

        int iang = 0;
        for( auto ang: ang_quad_ ) {

            // Configure the current worker for this angle
            cw.set_octant( iang / ang_quad_.ndir_oct() + 1 );

            real_t wgt = ang.weight * HPI; 
            real_t ox = ang.ox;
            real_t oy = ang.oy;
            real_t oz = ang.oz;

            // Configure the loop direction. Could template the below for speed
            // at some point.
            int sttx = 0;
            int stpx = nx_;
            int xdir = 1;
            if( ox < 0.0 ) {
                ox = -ox;
                sttx = nx_-1;
                stpx = -1;
                xdir = -1;
            }
            
            int stty = 0;
            int stpy = ny_;
            int ydir = 1;
            if( oy < 0.0 ) {
                oy = -oy;
                stty = ny_-1;
                stpy = -1;
                ydir = -1;
            }
            
            int sttz = 0;
            int stpz = nz_;
            int zdir = 1;
            if( oz < 0.0 ) {
                oz = -oz;
                sttz = nz_-1;
                stpz = -1;
                zdir = -1;
            }

            // initialize upwind condition
            x_flux = bc_in_.get_face( group, iang, Normal::X_NORM );
            y_flux = bc_in_.get_face( group, iang, Normal::Y_NORM );
            z_flux = bc_in_.get_face( group, iang, Normal::Z_NORM );

            cw.upwind_work( x_flux, y_flux, z_flux, ang, group);

            for( int iz=sttz; iz!=stpz; iz+=zdir ) {
                real_t tz = oz/hz_[iz];
                for( int iy=stty; iy!=stpy; iy+=ydir ) {
                    real_t ty = oy/hy_[iy];
                    for( int ix=sttx; ix!=stpx; ix+=xdir ) {
                        // Gross. really need an Sn mesh abstraction
                        real_t psi_lx = x_flux[ny_*iz + iy];
                        real_t psi_ly = y_flux[nx_*iz + ix];
                        real_t psi_lz = z_flux[nx_*iy + ix];

                        int i = iz*nx_*ny_ + iy*nx_ + ix;
                        real_t tx = ox/hx_[ix];
                        real_t psi = 2.0*(tx*psi_lx + 
                                           ty*psi_ly + 
                                           tz*psi_lz) + q_(i);
                        psi /= 2.0*(tx + ty + tz) + xstr_(i);

                        flux_1g_(i) += psi*wgt;

                        x_flux[ny_*iz + iy] = 2.0*psi - x_flux[ny_*iz + iy];
                        y_flux[nx_*iz + ix] = 2.0*psi - y_flux[nx_*iz + ix];
                        z_flux[nx_*iy + ix] = 2.0*psi - z_flux[nx_*iy + ix];

                        // Stash currents (or not, depending on the
                        // CurrentWorker template parameter)
                        cw.current_work( x_flux[ny_*iz + iy],
                                         y_flux[nx_*iz + ix],
                                         z_flux[nx_*iy + ix], i, ang, group );
                    }
                }

            }

            // store the downwind boundary condition
            bc_out_.set_face(0, iang, Normal::X_NORM, x_flux);
            bc_out_.set_face(0, iang, Normal::Y_NORM, y_flux);
            bc_out_.set_face(0, iang, Normal::Z_NORM, z_flux);
            iang++;
        }
        // Update the boundary condition
        this->update_boundary( group );

        return;
    }

    void SnSweeper::initialize() {
        flux_.fill(1.0);
        flux_old_.fill(1.0);
        bc_in_.initialize(1.0/FPI);
        
        return;
    }


    void SnSweeper::get_pin_flux_1g( int ig, VecF& flux ) const { 
        flux.resize( mesh_.n_pin() );

        for( size_t i=0; i<flux.size(); i++ ) {
            flux[i] = flux_(i, ig);
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
