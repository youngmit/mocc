#include "moc_sweeper_2d3d.hpp"

#include <algorithm>
#include <iostream>
#include <valarray>

using std::cout;
using std::cin;
using std::endl;

namespace mocc {
    MoCSweeper_2D3D::MoCSweeper_2D3D( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        MoCSweeper( input, mesh ),
        corrections_( nullptr )
    {
        
    };

    void MoCSweeper_2D3D::sweep1g_final( int group ) {
        flux_1g_.fill(0.0);
    
        ArrayX e_tau(rays_.max_segments(), 1);
        ArrayX psi1(rays_.max_segments()+1, 1);
        ArrayX psi2(rays_.max_segments()+1, 1);
        
        // Temporary storage for surface- and node-average fluxes
        ArrayF flux_surf_sum(mesh_.n_surf()*2);
        std::valarray<int> flux_surf_norm(mesh_.n_surf()*2);
        ArrayF flux_vol_sum(mesh_.n_pin()*2);
        ArrayF flux_vol_norm(mesh_.n_pin());
        ArrayF sigt_sum(mesh_.n_pin()*2);

    
        // Wipe out the existing currents
        coarse_data_->current.col( group ) = 0.0;
    
        int iplane = 0;
        for( auto &plane_rays: rays_ ) {
            int first_reg = mesh_.first_reg_plane(iplane);
            int iang = 0;
            // Angles
            for( auto &ang_rays: plane_rays ) {
                int iang1 = iang;
                int iang2 = ang_quad_.reverse(iang);
                Angle ang = ang_quad_[iang];

                // Zero out all of the flux sum arrays
                flux_surf_sum = 0.0;
                flux_surf_norm = 0.0;
                flux_vol_sum = 0.0;
                flux_vol_norm = 0.0;
                sigt_sum = 0.0;

                float_t stheta = sin(ang.theta);
                float_t rstheta = 1.0/stheta;
                float_t wt_v_st = ang.weight * rays_.spacing(iang) *
                    stheta * PI;

                float_t current_weights[3] = {
                    ang.weight * rays_.spacing(iang) * ang.ox,
                    ang.weight * rays_.spacing(iang) * ang.oy,
                    ang.weight * rays_.spacing(iang) * ang.oz
                };
    
                int iray = 0;
                for( auto &ray: ang_rays ) {
                    int bc1 = ray.bc(0);
                    int bc2 = ray.bc(1);
    
                    // Compute exponentials
                    for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                        int ireg = ray.seg_index(iseg) + first_reg;
                        e_tau(iseg) = 1.0 - exp( -xstr_(ireg) * 
                                ray.seg_len(iseg) * rstheta );
                    }
    
                    // Forward direction
                    // Initialize from bc
                    psi1(0) = boundary_[group][iplane][iang1][bc1];
    
                    // Propagate through core geometry
                    for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                        int ireg = ray.seg_index(iseg) + first_reg;
                        float_t psi_diff = (psi1(iseg) - qbar_(ireg)) * 
                        e_tau(iseg);
                        psi1(iseg+1) = psi1(iseg) - psi_diff;
                        flux_1g_(ireg) += psi_diff*wt_v_st;
                    }
                    // Store boundary condition
                    boundary_out_[iplane][iang1][bc2] = psi1( ray.nseg() );
    
                    // Backward direction
                    // Initialize from bc
                    psi2(ray.nseg()) =
                        boundary_[group][iplane][iang2][bc2];
    
                    // Propagate through core geometry
                    for( int iseg=ray.nseg()-1; iseg>=0; iseg-- ) {
                        int ireg = ray.seg_index(iseg) + first_reg;
                        float_t psi_diff = (psi2(iseg+1) - qbar_(ireg)) * 
                            e_tau(iseg);
                        psi2(iseg) = psi2(iseg+1) - psi_diff;
                        flux_1g_(ireg) += psi_diff*wt_v_st;
                    }
                    // Store boundary condition
                    boundary_out_[iplane][iang2][bc1] = psi2(0);

                    // Stash currents
                    size_t cell_fw = ray.cm_start_fw();
                    size_t cell_bw = ray.cm_start_bw();
                    size_t iseg_fw = 0;
                    size_t iseg_bw = ray.nseg();
                    auto begin = ray.cm_data().cbegin();
                    auto end = ray.cm_data().cend();
                    for( auto crd = begin; crd != end; ++crd ) {
                        size_t surf_fw = mesh_.coarse_surf( cell_fw, crd->fw );
                        size_t surf_bw = mesh_.coarse_surf( cell_bw, crd->bw );

                        Normal norm_fw = surface_to_normal( crd->fw );
                        Normal norm_bw = surface_to_normal( crd->bw );

                        coarse_data_->current( surf_fw, group ) += 
                            psi1(iseg_fw) * current_weights[(int)norm_fw];
                        coarse_data_->current( surf_bw, group ) -= 
                            psi2(iseg_bw) * current_weights[(int)norm_bw];

                        flux_surf_sum[surf_fw*2+0] += psi1(iseg_fw);
                        flux_surf_norm[surf_fw*2+0] += 1.0;
                        flux_surf_sum[surf_bw*2+1] += psi2(iseg_bw);
                        flux_surf_norm[surf_bw*2+1] += 1.0;

                        // Store foreward volumetric quantities
                        for( int i=0; i<crd->nseg_fw; i++ ) {
                            int ireg = ray.seg_index(iseg_fw) + first_reg;
                            float_t xstr = xstr_(ireg);
                            float_t t = ang.rsintheta * ray.seg_len(iseg_fw);
                            float_t fluxvol = t * qbar_(ireg) + e_tau(iseg_fw)*
                                (psi1(iseg_fw)-qbar_(ireg))/xstr;
                            flux_vol_sum[cell_fw*2+0] += fluxvol;
                            flux_vol_norm[cell_fw] += t;
                            sigt_sum[cell_fw*2+0] += xstr*fluxvol;
                            iseg_fw++;
                        }

                        // Store backward volumetric quantities
                        for( int i=0; i<crd->nseg_bw; i++ ) {
                            iseg_bw--;
                            int ireg = ray.seg_index(iseg_bw) + first_reg;
                            float_t xstr = xstr_(ireg);
                            float_t t = ang.rsintheta * ray.seg_len(iseg_bw);
                            float_t fluxvol = t * qbar_(ireg) + e_tau(iseg_bw)*
                                (psi2(iseg_bw)-qbar_(ireg))/xstr;
                            flux_vol_sum[cell_bw*2+1] += fluxvol;
                            sigt_sum[cell_bw*2+1] += xstr*fluxvol;
                        }
                        
                        if( crd != begin ) {
                            cell_fw = mesh_.coarse_neighbor( cell_fw, 
                                    (crd)->fw );
                            cell_bw = mesh_.coarse_neighbor( cell_bw, 
                                    (crd)->bw );
                        }

                    }
                    iray++;
                } // Rays

                // Normalize the flux and sigt values and calculate
                // correction factors for the current angle/energy
                for( size_t i=0; i<flux_surf_norm.size(); i++) {
                    flux_surf_sum[i] /= flux_surf_norm[i];
                }

                for( size_t i=0; i<flux_vol_norm.size(); i++) {
                    sigt_sum[2*i+0] /= flux_vol_sum[2*i+0];
                    sigt_sum[2*i+1] /= flux_vol_sum[2*i+1];
                    flux_vol_sum[2*i+0] /= flux_vol_norm[i];
                    flux_vol_sum[2*i+1] /= flux_vol_norm[i];
                }

                calculate_corrections( iang, group, flux_surf_sum, flux_vol_sum, 
                        sigt_sum );

                iang++;
            } // angles
            iplane++;
    
            // Scale the scalar flux by the volume and add back the source
            flux_1g_ = flux_1g_/(xstr_*vol_) + qbar_*FPI;
        } // planes
    
    
        this->update_boundary( group );
    
        return;
    }

    void MoCSweeper_2D3D::calculate_corrections( size_t ang, size_t group, 
            ArrayF flux_surf, ArrayF flux_node, ArrayF sigt ) {
        const int FW = 0;
        const int BW = 1;

        const int XL = 0;
        const int XR = 1;
        const int YL = 2;
        const int YR = 3;
        int iang1 = ang;
        int iang2 = ang_quad_.reverse(ang);
        float_t ox = ang_quad_[ang].ox;

        Surface surfs[2][4];
        // We know that all of our moc angles are positive in the y direction
        surfs[FW][YL] = Surface::SOUTH;
        surfs[FW][YR] = Surface::NORTH;
        surfs[BW][YL] = Surface::NORTH;
        surfs[BW][YR] = Surface::SOUTH;

        if( ox > 0.0 ) {
            surfs[FW][XL] = Surface::WEST;
            surfs[FW][XR] = Surface::EAST;
            surfs[BW][XL] = Surface::EAST;
            surfs[BW][XR] = Surface::WEST;
        } else {
            surfs[FW][XL] = Surface::EAST;
            surfs[FW][XR] = Surface::WEST;
            surfs[BW][XL] = Surface::WEST;
            surfs[BW][XR] = Surface::EAST;
        }

        for( size_t ic=0; ic<mesh_.n_pin(); ic++ ) {
            float xstr = (*sn_xs_mesh_)[ic].xsmactr()[group];

            // FW direction
            {
                float_t psi_xl = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[FW][XL])*2+0];
                float_t psi_xr = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[FW][XR])*2+0];
                float_t psi_yl = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[FW][YL])*2+0];
                float_t psi_yr = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[FW][YR])*2+0];
            
                float_t ax = flux_node[ic*2+0]/(psi_xl + psi_xr);
                float_t ay = flux_node[ic*2+0]/(psi_yl + psi_yr);

                float_t b = sigt[ic*2+0]/xstr;

                corrections_->alpha( ic, iang1, group, Normal::X_NORM ) = ax;
                corrections_->alpha( ic, iang1, group, Normal::Y_NORM ) = ay;
            
                corrections_->beta( ic, iang1, group ) = b;
            }

            // BW direction
            {
                float_t psi_xl = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[BW][XL])*2+1];
                float_t psi_xr = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[BW][XR])*2+1];
                float_t psi_yl = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[BW][YL])*2+1];
                float_t psi_yr = 
                    flux_surf[mesh_.coarse_surf(ic, surfs[BW][YR])*2+1];

                float_t ax = flux_node[ic*2+1]/(psi_xl + psi_xr);
                float_t ay = flux_node[ic*2+1]/(psi_yl + psi_yr);

                float_t b = sigt[ic*2+1]/xstr;

                corrections_->alpha( ic, iang2, group, Normal::X_NORM ) = ax;
                corrections_->alpha( ic, iang2, group, Normal::Y_NORM ) = ay;

                corrections_->beta( ic, iang1, group ) = b;
            }
        }

//        cout << "Surface fluxes: " << endl;
//        for( auto v: flux_surf ) {
//            cout << v << endl;
//        }
//
//        cout << "Cell fluxes: " << flux_node.size() << endl;
//        for( auto v: flux_node ) {
//            cout << v << endl;
//        }
//
//        cout << "sig t:" << endl;
//        for( auto v: sigt ) {
//            cout << v << endl;
//        }
        return;
    }

}
