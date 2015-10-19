#pragma once

#include <iostream>

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "coarse_data.hpp"
#include "core_mesh.hpp"
#include "eigen_interface.hpp"
#include "ray_data.hpp"
#include "transport_sweeper.hpp"
#include "xs_mesh.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class MoCSweeper: public TransportSweeper {
        struct BC {
            real_t fw;
            real_t bw;
        };
        typedef std::vector< // group
                std::vector< // plane
                std::vector< // angle
                std::vector< real_t > > > > BCSet_t;   // BCs
        typedef std::vector< // plane
                std::vector< // angle
                std::vector< real_t > > > BCSet_Out_t; // BCs
    public:
        MoCSweeper( const pugi::xml_node &input,
                    const CoreMesh& mesh );
        
        ~MoCSweeper() { }
        
        virtual void sweep(int group);

        void initialize();

        void get_pin_flux_1g( int group, VecF& flux ) const;

        void output( H5::CommonFG *node ) const;

        void homogenize( CoarseData &data ) const;

        /**
         * Return a copy of the sweeper's angular quadrature.
         */
        AngularQuadrature get_ang_quad() const {
            return ang_quad_;
        }

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return SP_XSMeshHomogenized_t( 
                    new XSMeshHomogenized( mesh_ ) );
        }
    protected:
        const CoreMesh& mesh_;

        AngularQuadrature ang_quad_;
        RayData rays_;
        
        // Boundary condition. ordered by energy, plane, angle, ray
        BCSet_t boundary_;
        BCSet_Out_t boundary_out_;

        // Array of one group transport cross sections
        ArrayF xstr_;

        // Temporary storage for 1-group scalar flux
        ArrayF flux_1g_;

        // One-group, isotropic source, scaled by transport cross section
        ArrayF qbar_;

        // Number of inner iterations per group sweep
        unsigned int n_inner_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;

        // Update the boundary conditions
        void update_boundary( int group );

        /**
         * \brief Perform an MoC sweep
         *
         * This method performs a single source iteration using MoC, sweeping
         * all angles and rays once for a given angle. The \c cw parameter
         * allows for auxiliary work to be done during the sweep, without
         * affecting runtime performance when not needed. Examples of this are
         * currents for CMFD coupling and correction factors for 2D3D/CDD
         * coupling. See \ref moc::Current and cmdo::CurrentCorrections
         * for examples of these.
         */
        template <typename CurrentWorker>
        void sweep1g( int group, CurrentWorker &cw ) {
            flux_1g_ = 0.0;
        
            ArrayF e_tau(rays_.max_segments());
            ArrayF psi1(rays_.max_segments()+1);
            ArrayF psi2(rays_.max_segments()+1);
        
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
        
                    // Set up the current worker for sweeping this angle
                    cw.set_angle( ang, rays_.spacing(iang) );
        
                    real_t stheta = sin(ang.theta);
                    real_t rstheta = 1.0/stheta;
                    real_t wt_v_st = ang.weight * rays_.spacing(iang) *
                        stheta * PI;
        
                    int iray = 0;
                    for( auto &ray: ang_rays ) {
                        int bc1 = ray.bc(0);
                        int bc2 = ray.bc(1);
        
                        // Compute exponentials
                        for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            e_tau[iseg] = 1.0 - exp( -xstr_[ireg] * 
                                    ray.seg_len(iseg) * rstheta );
                        }
        
                        // Forward direction
                        // Initialize from bc
                        psi1[0] = boundary_[group][iplane][iang1][bc1];
        
                        // Propagate through core geometry
                        for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            real_t psi_diff = (psi1[iseg] - qbar_[ireg]) * 
                            e_tau[iseg];
                            psi1[iseg+1] = psi1[iseg] - psi_diff;
                            flux_1g_[ireg] += psi_diff*wt_v_st;
                        }
                        // Store boundary condition
                        boundary_out_[iplane][iang1][bc2] = psi1[ ray.nseg() ];
        
                        // Backward direction
                        // Initialize from bc
                        psi2[ray.nseg()] =
                            boundary_[group][iplane][iang2][bc2];
        
                        // Propagate through core geometry
                        for( int iseg=ray.nseg()-1; iseg>=0; iseg-- ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            real_t psi_diff = (psi2[iseg+1] - qbar_[ireg]) * 
                                e_tau[iseg];
                            psi2[iseg] = psi2[iseg+1] - psi_diff;
                            flux_1g_[ireg] += psi_diff*wt_v_st;
                        }
                        // Store boundary condition
                        boundary_out_[iplane][iang2][bc1] = psi2[0];
        
                        // Stash currents
                        cw.post_ray( psi1, psi2, e_tau, ray, first_reg, group );
                        
                        iray++;
                    } // Rays
                    cw.post_angle( iang, group );
                    iang++;
                } // angles
                iplane++;
        
                // Scale the scalar flux by the volume and add back the source
                flux_1g_ = flux_1g_/(xstr_*vol_) + qbar_*FPI;
            } // planes
        
        
            this->update_boundary( group );
        
            return;
        }
    };
}
