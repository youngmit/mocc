#pragma once

#include <iostream>
#include <omp.h>

#include "pugixml.hpp"

#include "moc/ray_data.hpp"


#include "core/angular_quadrature.hpp"
#include "core/boundary_condition.hpp"
#include "core/coarse_data.hpp"
#include "core/core_mesh.hpp"
#include "core/eigen_interface.hpp"
#include "core/exponential.hpp"
#include "core/transport_sweeper.hpp"
#include "core/xs_mesh.hpp"
#include "core/xs_mesh_homogenized.hpp"

namespace mocc { namespace moc {
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

        void get_pin_flux_1g( int group, ArrayB1& flux ) const;

        real_t set_pin_flux_1g( int group, const ArrayB1 &pin_flux );

        void output( H5Node &node ) const;

        void homogenize( CoarseData &data ) const {
            throw EXCEPT("Not Implemented");
        }

        /**
         * \copybrief TransportSweeper::create_source()
         *
         * This mostly calls the base \ref TransportSweeper::create_source()
         * method, but also makes sure that the source is configured properly
         * for MoC.
         */
        virtual UP_Source_t create_source( const pugi::xml_node &input ) const {
            auto source = TransportSweeper::create_source( input );
            source->set_scale_transport(true);
            return source;
        }

        /**
         * Return a copy of the sweeper's angular quadrature.
         */
        AngularQuadrature get_ang_quad() const {
            return ang_quad_;
        }

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            auto xsm = SP_XSMeshHomogenized_t( new XSMeshHomogenized( mesh_ ) );
            xsm->set_flux(flux_);
            return xsm;
        }

        void check_balance( int group ) const;

    protected:
        // Data
        const CoreMesh& mesh_;

        RayData rays_;

        // Boundary condition. ordered by energy, plane, angle, ray
        BCSet_t boundary_;
        BCSet_Out_t boundary_out_;

        // Array of one group transport cross sections
        ArrayB1 xstr_;

        // Reference to a one-group slice of flux_. This should be
        // default-constructed, so that it only references data in flux_
        ArrayB1 flux_1g_;

        // Number of inner iterations per group sweep
        unsigned int n_inner_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;

        // Update the boundary conditions
        void update_boundary( int group );

        // Exponential table
        Exponential exp_;

        bool dump_rays_;

        /**
         * \brief Perform an MoC sweep
         *
         * This method performs a single source iteration using MoC, sweeping
         * all angles and rays once for a given angle. The \c cw parameter
         * allows for auxiliary work to be done during the sweep, without
         * affecting runtime performance when not needed. Examples of this are
         * currents for CMFD coupling and correction factors for 2D3D/CDD
         * coupling. See \ref moc::Current and \ref cmdo::CurrentCorrections
         * for examples of these.
         */
        template <typename CurrentWorker>
        void sweep1g( int group, CurrentWorker &cw ) {
            flux_1g_ = 0.0;

#pragma omp parallel default(shared)
            {
            ArrayB1 e_tau(rays_.max_segments());
            ArrayB1 psi1(rays_.max_segments()+1);
            ArrayB1 psi2(rays_.max_segments()+1);
            ArrayB1 t_flux(n_reg_);
            t_flux = 0.0;

            int iplane = 0;
            for( const auto plane_ray_id: mesh_.unique_planes() ) {
                cw.set_plane( iplane );
                const auto &plane_rays = rays_[plane_ray_id];
                int first_reg = mesh_.first_reg_plane(iplane);
                int iang = 0;
                // Angles
                for( const auto &ang_rays: plane_rays ) {
                    // Get the source for this angle
                    auto &qbar = source_->get_transport( iang );
                    int iang1 = iang;
                    int iang2 = ang_quad_.reverse(iang);
                    Angle ang = ang_quad_[iang];

                    // Set up the current worker for sweeping this angle
                    cw.set_angle( ang, rays_.spacing(iang) );

                    real_t stheta = sin(ang.theta);
                    real_t rstheta = 1.0/stheta;
                    real_t wt_v_st = ang.weight * rays_.spacing(iang) *
                        stheta * PI;

#pragma omp for schedule(static, 1)
                    for( int iray=0; iray<(int)ang_rays.size(); iray++ ) {
                        const auto& ray = ang_rays[iray];

                        int bc1 = ray.bc(0);
                        int bc2 = ray.bc(1);

                        // Compute exponentials
                        for( int iseg=0; iseg<ray.nseg(); iseg++ ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            e_tau(iseg) = 1.0 - exp_.exp( -xstr_(ireg) *
                                    ray.seg_len(iseg) * rstheta );
                        }

                        // Forward direction
                        // Initialize from bc
                        psi1(0) = boundary_[group][iplane][iang1][bc1];

                        // Propagate through core geometry
                        for( int iseg=0; iseg<ray.nseg(); iseg++ ) {
                            int ireg = ray.seg_index(iseg) + first_reg;
                            real_t psi_diff = (psi1(iseg) - qbar[ireg]) *
                                e_tau(iseg);
                            psi1(iseg+1) = psi1(iseg) - psi_diff;
                            t_flux(ireg) += psi_diff*wt_v_st;
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
                            real_t psi_diff = (psi2(iseg+1) - qbar[ireg]) *
                                e_tau(iseg);
                            psi2(iseg) = psi2(iseg+1) - psi_diff;
                            t_flux(ireg) += psi_diff*wt_v_st;
                        }
                        // Store boundary condition
                        boundary_out_[iplane][iang2][bc1] = psi2(0);

                        // Stash currents
                        cw.post_ray( psi1, psi2, e_tau, ray, first_reg, group );
                    } // Rays
                    cw.post_angle( iang, group );
                    iang++;
                } // angles
                iplane++;
            } // planes

#pragma omp barrier
#pragma omp critical
            {
                /// \todo make an array operation after refactoring out valarray
                for( int i=0; i<(int)n_reg_; i++ ) {
                    flux_1g_(i) += t_flux(i);
                }
            }
#pragma omp barrier
            // Scale the scalar flux by the volume and add back the source
#pragma omp single
            {
                // \todo this is not correct for angle-dependent sources!
                auto &qbar = source_->get_transport( 0 );
                for( int i=0; i<(int)n_reg_; i++)
                {
                    flux_1g_(i) = flux_1g_(i)/(xstr_(i)*vol_[i]) + qbar[i]*FPI;
                }
            } // OMP single

            cw.post_sweep( group );

            this->update_boundary( group );
            }

            return;
        }
    };
} }
