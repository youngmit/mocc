#pragma once

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "coarse_data.hpp"
#include "core_mesh.hpp"
#include "error.hpp"
#include "global_config.hpp"
#include "mesh.hpp"
#include "sn_boundary.hpp"
#include "sn_current_worker.hpp"
#include "sn_source.hpp"
#include "transport_sweeper.hpp"
#include "xs_mesh_homogenized.hpp"

namespace mocc {
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh );

        ~SnSweeper() { }

        virtual void sweep( int group ) = 0;

        void initialize();

        void get_pin_flux_1g( int ig, VecF& flux ) const;

        void output( H5::CommonFG *node ) const;

        // Override the create_source() method to make an SnSource instead of
        // the regular
        UP_Source_t create_source() const {
            Source *s = new SnSource( n_reg_, xs_mesh_.get(), this->cflux());
            UP_Source_t source( s );
            return source;
        }

        void homogenize( CoarseData &data ) const;

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return std::static_pointer_cast<XSMeshHomogenized>( xs_mesh_ );
        }

    protected:
        const CoreMesh &mesh_;

        // Update the boundary conditions 
        void update_boundary( int group );

        unsigned int n_inner_;
        AngularQuadrature ang_quad_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;
        
        // Temporary storage for 1-group scalar flux
        ArrayF flux_1g_;
        
        // Temporary storage of the current-group transport cross section
        ArrayF xstr_;

        // Single-group isotropic source, should include in-scatter
        ArrayF q_;

        // Incomming boundary condition
        SnBoundary bc_in_;

        // Outgoing boundary condition. Only difined for one group
        SnBoundary bc_out_;

        /**
         * \brief Generic Sn sweep procedure for orthogonal mesh. 
         *
         * This routine performs a single, one-group transport sweep with Sn. It
         * is templated on two parameters to tailor it to different differencing
         * schemes and current calculation requirements. To see examples of
         * template parameters, look at \ref sn::CellWorker and \ref
         * sn::Current.
         */
        template <typename CurrentWorker, typename CellWorker>
        void sweep_1g( int group, CellWorker &cell_worker ) {
            CurrentWorker cw( coarse_data_, &mesh_ );
            flux_1g_ = 0.0;

            int nx = mesh_.nx();
            int ny = mesh_.ny();
            int nz = mesh_.nz();
    
    		ArrayF x_flux(ny*nz);
    		ArrayF y_flux(nx*nz);
    		ArrayF z_flux(nx*ny);
    
            int iang = 0;
            for( auto ang: ang_quad_ ) {
                // Configure the current worker for this angle
                cw.set_octant( iang / ang_quad_.ndir_oct() + 1 );

                cell_worker.set_angle( iang, ang );
    
                real_t wgt = ang.weight * HPI; 
                real_t ox = ang.ox;
                real_t oy = ang.oy;
                real_t oz = ang.oz;
    
                // Configure the loop direction. Could template the below for speed
                // at some point.
                int sttx = 0;
                int stpx = nx;
                int xdir = 1;
                if( ox < 0.0 ) {
                    ox = -ox;
                    sttx = nx-1;
                    stpx = -1;
                    xdir = -1;
                }
                
                int stty = 0;
                int stpy = ny;
                int ydir = 1;
                if( oy < 0.0 ) {
                    oy = -oy;
                    stty = ny-1;
                    stpy = -1;
                    ydir = -1;
                }
                
                int sttz = 0;
                int stpz = nz;
                int zdir = 1;
                if( oz < 0.0 ) {
                    oz = -oz;
                    sttz = nz-1;
                    stpz = -1;
                    zdir = -1;
                }
    
                // initialize upwind condition
                x_flux = bc_in_.get_face( group, iang, Normal::X_NORM );
                y_flux = bc_in_.get_face( group, iang, Normal::Y_NORM );
                z_flux = bc_in_.get_face( group, iang, Normal::Z_NORM );
    
                cw.upwind_work( x_flux, y_flux, z_flux, ang, group);
    
                for( int iz=sttz; iz!=stpz; iz+=zdir ) {
                    cell_worker.set_z(iz);
                    for( int iy=stty; iy!=stpy; iy+=ydir ) {
                        cell_worker.set_y(iy);
                        for( int ix=sttx; ix!=stpx; ix+=xdir ) {
                            // Gross. really need an Sn mesh abstraction
                            real_t psi_x = x_flux[ny*iz + iy];
                            real_t psi_y = y_flux[nx*iz + ix];
                            real_t psi_z = z_flux[nx*iy + ix];
    
                            int i = mesh_.coarse_cell( Position( ix, iy, iz ) );
                            
                            real_t psi = cell_worker.evaluate( psi_x, psi_y,
                                    psi_z, q_[i], xstr_[i], i );
    
                            x_flux[ny*iz + iy] = psi_x;
                            y_flux[nx*iz + ix] = psi_y;
                            z_flux[nx*iy + ix] = psi_z;
    
                            flux_1g_[i] += psi*wgt;

                            // Stash currents (or not, depending on the
                            // CurrentWorker template parameter)
                            cw.current_work( x_flux[ny*iz + iy],
                                             y_flux[nx*iz + ix],
                                             z_flux[nx*iy + ix], i, ang, group );
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
    };

    typedef std::unique_ptr<SnSweeper> UP_SnSweeper_t;
}
