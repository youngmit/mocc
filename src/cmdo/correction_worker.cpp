#include "correction_worker.hpp"

namespace mocc {
    namespace cmdo {
        void CurrentCorrections::calculate_corrections( size_t ang, 
                size_t group ) 
        {
            const int FW = 0;
            const int BW = 1;
            
            const int XL = 0;
            const int XR = 1;
            const int YL = 2;
            const int YR = 3;
            int iang1 = ang;
            int iang2 = ang_quad_.reverse(ang);
            real_t ox = ang_quad_[ang].ox;
            
            Surface surfs[2][4];
            // We know that all of our moc angles are positive in the y
            // direction
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
            
            // See the Surface Normalization page
            real_t area[2] = 
            {
                std::abs( rays_.spacing( ang )/cos(ang_quad_[ang].alpha) ),
                std::abs( rays_.spacing( ang )/sin(ang_quad_[ang].alpha) )
            };
            
            for( size_t ic=0; ic<mesh_->n_pin(); ic++ ) {
                auto pos = mesh_->coarse_position(ic);
            
                real_t area_x = area[0]/mesh_->pin_dx()[pos.x];
                real_t area_y = area[1]/mesh_->pin_dy()[pos.y];

                real_t xstr = sn_xs_mesh_[ic].xsmactr()[group];
            
                // FW direction
                {
                    real_t psi_xl = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[FW][XL])*2+0]*area_x;
                    real_t psi_xr = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[FW][XR])*2+0]*area_x;
                    real_t psi_yl = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[FW][YL])*2+0]*area_y;
                    real_t psi_yr = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[FW][YR])*2+0]*area_y;
                
            
                    real_t ax = vol_sum_[ic*2+0]/(psi_xl + psi_xr);
                    real_t ay = vol_sum_[ic*2+0]/(psi_yl + psi_yr);
            
                    real_t b = sigt_sum_[ic*2+0]/xstr;
            
                    corrections_->alpha( ic, iang1, group, Normal::X_NORM ) = ax;
                    corrections_->alpha( ic, iang1, group, Normal::Y_NORM ) = ay;
                
                    corrections_->beta( ic, iang1, group ) = b;
                }
            
                // BW direction
                {
                    real_t psi_xl = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[BW][XL])*2+1]*area_x;
                    real_t psi_xr = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[BW][XR])*2+1]*area_x;
                    real_t psi_yl = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[BW][YL])*2+1]*area_y;
                    real_t psi_yr = 
                        surf_sum_[mesh_->coarse_surf(ic, surfs[BW][YR])*2+1]*area_y;
            
                    real_t ax = vol_sum_[ic*2+1]/(psi_xl + psi_xr);
                    real_t ay = vol_sum_[ic*2+1]/(psi_yl + psi_yr);
            
                    real_t b = sigt_sum_[ic*2+1]/xstr;
            
                    corrections_->alpha( ic, iang2, group, Normal::X_NORM ) = ax;
                    corrections_->alpha( ic, iang2, group, Normal::Y_NORM ) = ay;
            
                    corrections_->beta( ic, iang2, group ) = b;
                }
            }

            return;
        }
    }
}
