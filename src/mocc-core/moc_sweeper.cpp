#include "moc_sweeper.hpp"

#include "error.hpp"

using std::endl;
using std::cout;

namespace mocc {
    
    MoCSweeper::MoCSweeper( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        TransportSweeper( mesh ),
        ang_quad_( input.child("ang_quad") ),
        rays_( input.child("rays"), ang_quad_, mesh ),
        xstr_( n_reg_, 1 ),
        flux_1g_( n_reg_, 1 ),
        qbar_( n_reg_, 1 ),
        bc_type_( mesh_.boundary() )
    {   

        // Make sure we have input from the XML
        if( input.empty() ) {
            Error("No input specified to initialize MoC sweeper.");
        }

        // Parse the number of inner iterations
        int int_in = input.attribute("n_inner").as_int(-1);
        if(int_in < 0) {
            Error("Invalid number of inner iterations specified (n_inner).");
        }
        n_inner_ = int_in;

        // Set up the array of volumes (surface area)
        int ireg = 0;
        for( auto pin=mesh_.begin_pin(); pin!=mesh_.end_pin(); ++pin ) {
            const VecF& pin_vol = (*pin)->vol();
            for( auto &v: pin_vol ) {
                vol_(ireg) = v;
                ireg++;
            }
        }

        // allocate space to store the boundary conditions
        boundary_.resize( ng_ );
        for( auto &group_rays: boundary_ ) {
            group_rays.resize( mesh_.nz() );
            for( auto &angle_rays: group_rays ) {
                // We actually allocate BCs for all 4 octants to make things a
                // little simpler.
                angle_rays.resize( ang_quad_.ndir_oct()*4 );
                int iang = 0;
                for( auto ang_it=ang_quad_.octant(1); 
                        ang_it!=ang_quad_.octant(5); ang_it++ ) {
                    angle_rays[iang].resize(rays_.n_rays(iang));
                    iang++;
                }
            }
        }
        boundary_out_.resize( mesh_.nz() );
        for( auto &angle_rays: boundary_out_ ) {
            // We actually allocate BCs for all 4 octants to make things a
            // little simpler.
            angle_rays.resize( ang_quad_.ndir_oct()*4 );
            int iang = 0;
            for( auto ang_it=ang_quad_.octant(1); 
                    ang_it!=ang_quad_.octant(5); ang_it++ ) {
                angle_rays[iang].resize(rays_.n_rays(iang));
                iang++;
            }
        }

        return;
    }

    void MoCSweeper::sweep( int group ) {
        // set up the xstr_ array
        for( auto &xsr: xs_mesh_ ) {
            float_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_(ireg) = xstr;
            }
        }

        flux_1g_ = flux_.col( group );
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            // update the self-scattering source
            source_->self_scatter( group, flux_1g_, qbar_ );
            this->sweep1g( group );
        }

        flux_.col( group ) = flux_1g_;

        return;
    }

    void MoCSweeper::sweep1g( int group ) {
        flux_1g_.fill(0.0);

        ArrayX e_tau(rays_.max_segments(), 1);

        int iplane = 0;
        for( auto &plane_rays: rays_ ) {
            int first_reg = mesh_.first_reg_plane(iplane);
            int iang = 0;
            for( auto &ang_rays: plane_rays ) {
                int iang1 = iang;
                int iang2 = ang_quad_.reverse(iang);
                float_t stheta = sin(ang_quad_[iang].theta);
                float_t rstheta = 1.0/stheta;
                float_t wt_v_st = ang_quad_[iang].weight * rays_.spacing(iang) *
                    stheta * PI;

                int iray = 0;
                for( auto &ray: ang_rays ) {
                    int bc1 = ray.bc(0);
                    int bc2 = ray.bc(1);

                    // Compute exponentials
                    for( int iseg=0; iseg<ray.nseg(); iseg++ ) {
                        int ireg = ray.seg_index(iseg) + first_reg;
                        e_tau(iseg) = 1.0 - exp( -xstr_(ireg) * 
                                ray.seg_len(iseg) * rstheta );
                    }


                    // Forward direction
                    {
                        // Initialize from bc
                        float_t psi = 
                            boundary_[group][iplane][iang1][bc1];

                        // Propagate through core geometry
                        for( int iseg=0; iseg<ray.nseg(); iseg++ ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            float_t psi_diff = (psi - qbar_(ireg)) * e_tau(iseg);
                            psi -= psi_diff;
                            flux_1g_(ireg) += psi_diff*wt_v_st;
                        }
                        // Store boundary condition
                        boundary_out_[iplane][iang1][bc2] = psi;
                    }
                    
                    // Backward direction
                    {
                        // Initialize from bc
                        float_t psi =
                            boundary_[group][iplane][iang2][bc2];

                        // Propagate through core geometry
                        for( int iseg=ray.nseg()-1; iseg>=0; iseg-- ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            float_t psi_diff = (psi - qbar_(ireg)) * e_tau(iseg);
                            psi -= psi_diff;
                            flux_1g_(ireg) += psi_diff*wt_v_st;
                        }
                        // Store boundary condition
                        boundary_out_[iplane][iang2][bc1] = psi;
                    }
                    iray++;
                } // Rays
                iang++;
            } // angles
            iplane++;

            // Scale the scalar flux by the volume and add back the source
            flux_1g_ = flux_1g_/(xstr_*vol_) + qbar_*FPI;
        } // planes


        this->update_boundary( group );

        return;
    }

    void MoCSweeper::update_boundary( int group ) {
        int iplane=0;
        for( auto &plane_bcs: boundary_[group] ) {
            for( unsigned int iang=0; iang<plane_bcs.size(); iang++ ) {
                int nx = rays_.nx(iang);
                int ny = rays_.ny(iang);

                // Determine based on the angle quadrant which surfaces are
                // upwind and need to be updated for the given angle
                Surface upwind[2];
                if(ang_quad_[iang].ox > 0.0) {
                    upwind[0] = WEST;
                } else {
                    upwind[0] = EAST;
                }
                if(ang_quad_[iang].oy > 0.0) {
                    upwind[1] = SOUTH;
                } else {
                    upwind[1] = NORTH;
                }

                // Treat the two surfaces that were determined above
                if( bc_type_[upwind[0]] == REFLECT ) {
                    int ang_ref = ang_quad_.reflect(iang, upwind[0]);
                    for( int ibc=0; ibc<ny; ibc++ ) {
                        plane_bcs[iang][ibc] =
                            boundary_out_[iplane][ang_ref][ibc];
                    }
                } else {
                    for( int ibc=0; ibc<ny; ibc++ ) {
                        plane_bcs[iang][ibc] = 0.0;
                    }
                }
                
                if( bc_type_[upwind[1]] == REFLECT ) {
                    int ang_ref = ang_quad_.reflect(iang, upwind[1]);
                    for( int ibc=ny; ibc<(nx+ny); ibc++ ) {
                        plane_bcs[iang][ibc] =
                            boundary_out_[iplane][ang_ref][ibc];
                    }
                } else {
                    for( int ibc=ny; ibc<(nx+ny); ibc++ ) {
                        plane_bcs[iang][ibc] = 0.0;
                    }
                }
            } // angle loop
            iplane++;
        } // plane loop
    }

    void MoCSweeper::initialize() {
        // There are better ways to do this, but for now, just start with 1.0
        flux_.fill(1.0);

        // Walk through the boundary conditions and initialize them the 1/4pi
        float_t val = 1.0/FPI;
        for( auto &group_rays: boundary_ ) {
            for( auto &plane_rays: group_rays ) {
                for( auto &angle_rays: plane_rays ) {
                    for( auto &ray: angle_rays ) {
                        ray = val;
                    }
                }
            }
        }
        return;
    }

    void MoCSweeper::calc_fission_source( float_t k, 
            ArrayX& fission_source ) const {

        float_t rkeff = 1.0/k;
        fission_source.fill(0.0);
        for( auto &xsr: xs_mesh_ ) {
            const auto& xsnf = xsr.xsmacnf();
            for(unsigned int ig=0; ig<xs_mesh_.n_grp(); ig++ ) {
                for( auto &ireg: xsr.reg() ) {
                    fission_source(ireg) += rkeff*xsnf[ig]*flux_(ireg, ig);
                }
            }
        }
        return;
    }

    void MoCSweeper::get_pin_flux( int group, VecF& flux) const {
        flux.resize( mesh_.nx()*mesh_.ny()*mesh_.nz() );
        for( auto &f: flux ) {
            f = 0.0;
        }

        int ireg = 0;
        int ipin = 0;
        for( auto &pin: mesh_ ) {
            Position pos = mesh_.pin_position(ipin);
            int i = mesh_.index_lex(pos);

            for( int ir=0; ir<pin->n_reg(); ir++) {
                flux[i] += flux_(ireg, group)*vol_(ireg);
                ireg++;
            }

            ipin++;
        }

        return;
    }
}
