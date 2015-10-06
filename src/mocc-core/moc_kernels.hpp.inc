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
                for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
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
                    for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
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

void MoCSweeper::sweep1g_final( int group ) {
    flux_1g_.fill(0.0);

    ArrayX e_tau(rays_.max_segments(), 1);
    ArrayX psi1(rays_.max_segments()+1, 1);
    ArrayX psi2(rays_.max_segments()+1, 1);

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



                size_t cell_fw = ray.cm_cell_fw();
                size_t cell_bw = ray.cm_cell_bw();
                size_t surf_fw = ray.cm_surf_fw();
                size_t surf_bw = ray.cm_surf_bw();
                size_t iseg_fw = 0;
                size_t iseg_bw = ray.nseg();

                Normal norm_fw = mesh_.surface_normal( surf_fw );
                Normal norm_bw = mesh_.surface_normal( surf_bw );
                coarse_data_->current( surf_fw, group ) += 
                    psi1(iseg_fw) * current_weights[(int)norm_fw];
                coarse_data_->current( surf_bw, group ) -= 
                    psi2(iseg_bw) * current_weights[(int)norm_bw];
                auto begin = ray.cm_data().cbegin();
                auto end = ray.cm_data().cend();
                for( auto crd = begin; crd != end; ++crd ) {
                    // Hopefully branch prediction saves me here.
                    if( crd->fw != Surface::INVALID ) {
                        iseg_fw += crd->nseg_fw;
                        norm_fw = surface_to_normal( crd->fw );
                        surf_fw = mesh_.coarse_surf( cell_fw, crd->fw );
                        coarse_data_->current( surf_fw, group ) += 
                            psi1(iseg_fw) * current_weights[(int)norm_fw];
                    }

                    if( crd->bw != Surface::INVALID ) {
                        iseg_bw -= crd->nseg_bw;
                        norm_bw = surface_to_normal( crd->bw );
                        surf_bw = mesh_.coarse_surf( cell_bw, crd->bw );
                        coarse_data_->current( surf_bw, group ) -= 
                            psi2(iseg_bw) * current_weights[(int)norm_bw];
                    }

                    cell_fw = mesh_.coarse_neighbor( cell_fw, (crd)->fw );
                    cell_bw = mesh_.coarse_neighbor( cell_bw, (crd)->bw );
                }
                iray++;
            } // Rays
            iang++;
        } // angles
        iplane++;

//cout << endl;
//cout << coarse_data_->current.col(group) << endl;
//cout << endl;

//cin.ignore();

        // Scale the scalar flux by the volume and add back the source
        flux_1g_ = flux_1g_/(xstr_*vol_) + qbar_*FPI;
    } // planes


    this->update_boundary( group );

    return;
}
