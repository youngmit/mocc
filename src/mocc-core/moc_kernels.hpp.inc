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

/// \todo actually implement current accumulation
void MoCSweeper::sweep1g_current( int group ) {
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
