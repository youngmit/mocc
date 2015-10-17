#include "sn_current_worker.hpp"
#include "sn_sweeper_cdd.hpp"

using std::cout;
using std::cin;
using std::endl;

namespace mocc {
    SnSweeper_CDD::SnSweeper_CDD( const pugi::xml_node &input, 
            const CoreMesh &mesh):
        SnSweeper( input, mesh ),
        corrections_( nullptr )
    {
        if( input.child("data") ) {
            if( input.child("data").attribute("type") ) {
                std::string data_type = 
                    input.child("data").attribute("type").value();
                if( data_type == "default" ) {
                    cout << "Generating default values for correction factors."
                        << endl;
                    my_corrections_.reset( new CorrectionData( n_reg_, 
                                ang_quad_.ndir(), n_group_) );
                    corrections_ = my_corrections_.get();
                }
            }
        }
    }

    void SnSweeper_CDD::sweep( int group ) {
        // Store the transport cross section somewhere useful
        for( auto &xsr: *xs_mesh_ ) {
            real_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_[ireg] = xstr;
            }
        }
        
        flux_1g_ = flux_[std::slice(group*n_reg_, n_reg_, 1)];

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {

            // Set the source (add upscatter and divide by 4PI)
            source_->self_scatter( group, flux_1g_, q_ );

            if( inner == n_inner_-1 && coarse_data_ ) {
cout << "sn source" << endl;
for( auto q: q_ ) {
    cout << q << endl;
}
                this->sweep_std<sn::Current>( group );
            } else {
                this->sweep_std<sn::NoCurrent>( group );
            }
        }
        flux_[std::slice(group*n_reg_, n_reg_, 1)] = flux_1g_;

        return;
    }

    // Perform a single sweep as fast as possible with the CDD equations. Don't
    // collect currents or anything.
    template <typename CurrentWorker>
    void SnSweeper_CDD::sweep_std( int group ) {
        CurrentWorker cw( coarse_data_, &mesh_ );
        assert( corrections_ );
        flux_1g_ = 0.0;

		ArrayF x_flux(ny_*nz_);
		ArrayF y_flux(nx_*nz_);
		ArrayF z_flux(nx_*ny_);

        int nang_half = ang_quad_.ndir()/2;

        int iang = 0;
        for( auto ang: ang_quad_ ) {
            int iang_a = iang % nang_half;
            real_t wgt = ang.weight * HPI; 
            real_t ox = ang.ox;
            real_t oy = ang.oy;
            real_t oz = ang.oz;

            // Configure the loop direction. Could template this for speed at
            // some point.
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
                real_t tz = oz/hz_[iz];
                for( int iy=stty; iy!=stpy; iy+=ydir ) {
                    real_t ty = oy/hy_[iy];
                    for( int ix=sttx; ix!=stpx; ix+=xdir ) {
                        // Gross. really need an Sn mesh abstraction
                        int i = iz*nx_*ny_ + iy*nx_ + ix;
                        real_t tx = ox/hx_[ix];

                        real_t psi_lx = x_flux[ny_*iz + iy];
						real_t psi_ly = y_flux[nx_*iz + ix];
						real_t psi_lz = z_flux[nx_*iy + ix];

                        real_t ax = corrections_->alpha( i, iang_a, group, 
                                Normal::X_NORM);
                        real_t ay = corrections_->alpha( i, iang_a, group, 
                                Normal::Y_NORM);
                        real_t b = corrections_->beta( i, iang_a, group );

                        real_t gx = ax*b;
                        real_t gy = ay*b;

                        real_t psi = q_[i] + 
                            2.0*(tx*psi_lx + ty*psi_ly + tz*psi_lz );
                        psi /= tx/gx + ty/gy + 2.0*tz + xstr_[i];

                        flux_1g_[i] += psi*wgt;

                        x_flux[ny_*iz + iy] = (psi - gx*psi_lx) / gx;
						y_flux[nx_*iz + ix] = (psi - gy*psi_ly) / gy;
						z_flux[nx_*iy + ix] = 2.0*psi - psi_lz;

                        // Do current worker (maybe)
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
        } // angle loop

        // Update the boundary condition
        this->update_boundary( group );

        return;
    }
}
