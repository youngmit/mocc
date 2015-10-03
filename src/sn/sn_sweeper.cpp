#include "sn_sweeper.hpp"

#include <iomanip>

#include "error.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
    SnSweeper::SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh ):
        core_mesh_( mesh ),
        ang_quad_( input.child("ang_quad") ),
        bc_type_( mesh.boundary() ),
        xstr_( mesh.n_pin(), 1 ), 
        q_( mesh.n_pin(), 1 )
    {
        // Set up all of the stuff that would normally be done by the
        // TransportSweeper constructor. There is probably a better and more
        // maintainable way to do this; will revisit.
        xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh) );
        n_reg_ = mesh.n_pin();
        ng_ = xs_mesh_->n_group();
        flux_.resize( n_reg_, ng_ );
        flux_old_.resize( n_reg_, ng_ );
        vol_.resize( n_reg_, 1 );

        flux_1g_.resize( n_reg_, 1 );

        // Set the mesh volumes. Same as the pin volumes
        int ipin = 0;
        for( auto &pin: core_mesh_ ) {
            int i = core_mesh_.index_lex( core_mesh_.pin_position(ipin) );
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
        const VecF& core_x = core_mesh_.pin_hx();
        for( unsigned int i=0; i<core_x.size()-1; i++ ) {
            hx_.push_back(core_x[i+1]-core_x[i]);
        }
        const VecF& core_y = core_mesh_.pin_hy();
        for( unsigned int i=0; i<core_y.size()-1; i++ ) {
            hy_.push_back(core_y[i+1]-core_y[i]);
        }
        hz_ = core_mesh_.hz();

        nx_ = hx_.size();
        ny_ = hy_.size();
        nz_ = hz_.size();

        bc_in_ = SnBoundary( ng_, ang_quad_.ndir(), nx_, ny_, nz_ );
        bc_out_ = SnBoundary( 1, ang_quad_.ndir(), nx_, ny_, nz_ );

        return;
    }

    void SnSweeper::homogenize( CoarseData &data ) const {
        return;
    }

    void SnSweeper::sweep( int group ) {
        // Store the transport cross section somewhere useful
        for( auto &xsr: *xs_mesh_ ) {
            float_t xstr = xsr.xsmactr()[group];
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
                this->sweep_final( group );
            } else {
                this->sweep_std( group );
            }
        }
        flux_.col( group ) = flux_1g_;

        return;
    }

    void SnSweeper::sweep_std( int group ) {
        flux_1g_.fill(0.0);

		ArrayF x_flux(ny_*nz_);
		ArrayF y_flux(nx_*nz_);
		ArrayF z_flux(nx_*ny_);

        int iang = 0;
        for( auto ang: ang_quad_ ) {
            float_t wgt = ang.weight * HPI; 
            float_t ox = ang.ox;
            float_t oy = ang.oy;
            float_t oz = ang.oz;

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
            x_flux = bc_in_.get_face( group, iang, Normal::X_NORM);
            y_flux = bc_in_.get_face( group, iang, Normal::Y_NORM);
            z_flux = bc_in_.get_face( group, iang, Normal::Z_NORM);

            for( int iz=sttz; iz!=stpz; iz+=zdir ) {
                float_t tz = oz/hz_[iz];
                for( int iy=stty; iy!=stpy; iy+=ydir ) {
                    float_t ty = oy/hy_[iy];
                    for( int ix=sttx; ix!=stpx; ix+=xdir ) {
                        // Gross. really need an Sn mesh abstraction
                        float_t psi_lx = x_flux[ny_*iz + iy];
                        float_t psi_ly = y_flux[nx_*iz + ix];
                        float_t psi_lz = z_flux[nx_*iy + ix];

                        int i = iz*nx_*ny_ + iy*nx_ + ix;
                        float_t tx = ox/hx_[ix];
                        float_t psi = 2.0*(tx*psi_lx + 
                                           ty*psi_ly + 
                                           tz*psi_lz) + q_(i);
                        psi /= 2.0*(tx + ty + tz) + xstr_(i);

                        flux_1g_(i) += psi*wgt;

                        x_flux[ny_*iz + iy] = 2.0*psi - x_flux[ny_*iz + iy];
                        y_flux[nx_*iz + ix] = 2.0*psi - y_flux[nx_*iz + ix];
                        z_flux[nx_*iy + ix] = 2.0*psi - z_flux[nx_*iy + ix];
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

    /// \todo implement current calculations here
    void SnSweeper::sweep_final( int group ) {
        this->sweep_std( group );
        return;
    }

    void SnSweeper::initialize() {
        flux_.fill(1.0);
        flux_old_.fill(1.0);
        bc_in_.initialize(0.0);
        
        return;
    }


    void SnSweeper::get_pin_flux( int ig, VecF& flux ) const { 
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
                bc_in_.set_face(group, iang, norm, 
                        bc_out_.get_face(0, iang_refl, norm ));
            } else {
                bc_in_.zero_face(group, iang, norm);
            }
            // Y-normal surface
            norm = Normal::Y_NORM;
            surf = (ang.oy > 0) ? Surface::SOUTH : Surface::NORTH;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_type_[(int)surf] == Boundary::REFLECT ) {
                bc_in_.set_face(group, iang, norm, 
                        bc_out_.get_face(0, iang_refl, norm ));
            } else {
                bc_in_.zero_face(group, iang, norm);
            }
            // Z-normal surface
            norm = Normal::Z_NORM;
            surf = (ang.oz > 0) ? Surface::BOTTOM : Surface::TOP;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_type_[(int)surf] == Boundary::REFLECT ) {
                bc_in_.set_face(group, iang, norm, 
                        bc_out_.get_face(0, iang_refl, norm ));
            } else {
                bc_in_.zero_face(group, iang, norm);
            }

            iang++;
        }
        return;
    }

    void SnSweeper::output( H5File& file ) const {
        VecI dims;
        dims.push_back(nz_);
        dims.push_back(ny_);
        dims.push_back(nx_);
        
        // Make a group in the file to store the flux
        file.mkdir("/flux");
        
        // Provide energy group upper bounds
        file.write("/eubounds", xs_mesh_->eubounds(), VecI(1, ng_));
        file.write("/ng", ng_);
        
        for( unsigned int ig=0; ig<ng_; ig++ ) {
            VecF flux;
            
            for( unsigned int ireg=0; ireg<n_reg_; ireg++ ) {
                flux.push_back( flux_( ireg, ig ) );
            }
        
            std::stringstream setname;
            setname << "/flux/" << std::setfill('0') << std::setw(3) << ig+1;
        
            file.write(setname.str(), flux, dims);
        }
        return;
    }
}
