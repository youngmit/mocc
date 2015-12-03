#pragma once

#include "pugixml.hpp"

#include "mocc-core/angular_quadrature.hpp"
#include "mocc-core/blitz_typedefs.hpp"
#include "mocc-core/coarse_data.hpp"
#include "mocc-core/core_mesh.hpp"
#include "mocc-core/error.hpp"
#include "mocc-core/files.hpp"
#include "mocc-core/global_config.hpp"
#include "mocc-core/mesh.hpp"
#include "mocc-core/transport_sweeper.hpp"
#include "mocc-core/utils.hpp"
#include "mocc-core/xs_mesh_homogenized.hpp"

#include "sn/sn_boundary.hpp"
#include "sn/sn_current_worker.hpp"
#include "sn/sn_source.hpp"

namespace mocc { namespace sn {
    template <class Worker>
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node& input, const CoreMesh& mesh ):
                mesh_( mesh ),
                ang_quad_( input.child("ang_quad") ),
                bc_type_( mesh.boundary() ),
                flux_1g_( mesh_.n_pin() ),
                xstr_( mesh.n_pin() ),
                q_( mesh_.n_pin() ),
                bc_in_( mesh.mat_lib().n_group(), ang_quad_, mesh_ ),
                bc_out_( 1, ang_quad_, mesh_ ),
                cell_worker_( mesh_, ang_quad_ )
        {
            LogFile << "Constructing a base Sn sweeper" << std::endl;

            // Set up all of the stuff that would normally be done by the
            // TransportSweeper constructor. There is probably a better and more
            // maintainable way to do this; will revisit.
            core_mesh_ = &mesh;
            xs_mesh_ = SP_XSMesh_t( new XSMeshHomogenized(mesh) );
            n_reg_ = mesh.n_pin();
            n_group_ = xs_mesh_->n_group();
            flux_.resize( n_reg_, n_group_ );
            flux_old_.resize( n_reg_, n_group_ );
            vol_.resize( n_reg_ );

            // Set the mesh volumes. Same as the pin volumes
            int ipin = 0;
            for( auto &pin: mesh_ ) {
                int i = mesh_.index_lex( mesh_.pin_position(ipin) );
                vol_[i] = pin->vol();
                ipin++;
            }

            // Make sure we have input from the XML
            if( input.empty() ) {
                throw EXCEPT("No input specified to initialize Sn sweeper.");
            }

            // Parse the number of inner iterations
            int int_in = input.attribute("n_inner").as_int(-1);
            if( int_in < 0 ) {
                throw EXCEPT("Invalid number of inner iterations specified "
                        "(n_inner).");
            }
            n_inner_ = int_in;

            return;
        }

        ~SnSweeper() { }
        
        void initialize() {
            flux_ = 1.0;
            flux_old_ = 1.0;
            bc_in_.initialize(1.0/FPI);

            return;
        }

        void get_pin_flux_1g( int ig, ArrayB1& flux ) const {
            assert( flux.size() == n_reg_ );

            flux = flux_(blitz::Range::all(), ig);

            return;
        }

        /**
         * Just copy the flux across, since no homogenization is necessary.
         */
        real_t set_pin_flux_1g( int group, const ArrayB1 &pin_flux ) {
            assert( pin_flux.size() == n_reg_ );

            real_t resid = 0.0;
            size_t i = 0;
            for( auto &v: pin_flux ) {
                real_t e = flux_(1, group) - v;
                resid += e*e;
                flux_((int)i, (int)group) = v;
                i++;
            }
            return std::sqrt(resid);
        }

        /**
         * \brief Re-assign the angular quadrature.
         */
        void set_ang_quad( AngularQuadrature ang_quad ) {
            ang_quad_ = ang_quad;
            return;
        }

        Worker& get_worker() {
            return cell_worker_;
        }

        
        // Override the create_source() method to make an SnSource instead of
        // the regular
        UP_Source_t create_source() const {
            Source *s = new SnSource( n_reg_, xs_mesh_.get(), this->flux());
            UP_Source_t source( s );
            return source;
        }

        


        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return std::static_pointer_cast<XSMeshHomogenized>( xs_mesh_ );
        }

        
        void sweep( int group ) {
            /// \todo add an is_read() method to the worker, and make sure that
            /// it's ready before continuing.

            // Store the transport cross section somewhere useful
            for( auto &xsr: *xs_mesh_ ) {
                real_t xstr = xsr.xsmactr()[group];
                for( auto &ireg: xsr.reg() ) {
                    xstr_[ireg] = xstr;
                }
            }
    
            flux_1g_ = flux_(blitz::Range::all(), group);
    
            // Perform inner iterations
            for( size_t inner=0; inner<n_inner_; inner++ ) {
                // Set the source (add upscatter and divide by 4PI)
                source_->self_scatter( group, flux_1g_, q_ );
                if( inner == n_inner_-1 && coarse_data_ ) {
                    // Wipe out the existing currents
                    coarse_data_->current( blitz::Range::all(), group ) = 0.0;
                    this->sweep_1g<sn::Current>( group );
std::cout << coarse_data_->current(blitz::Range::all(), group) << std::endl;
                } else {
                    this->sweep_1g<sn::NoCurrent>( group );
                }
            }
            flux_(blitz::Range::all(), group) = flux_1g_;
    
            return;
        }

        void output( H5::CommonFG *node ) const {
            auto dims = mesh_.dimensions();
            std::reverse( dims.begin(), dims.end() );

            // Make a group in the file to store the flux
            node->createGroup("flux");

            ArrayB2 flux = this->get_pin_flux();
            Normalize( flux.begin(), flux.end() );

            for( unsigned int ig=0; ig<n_group_; ig++ ) {
                std::stringstream setname;
                setname << "flux/" << std::setfill('0') << std::setw(3) << ig+1;

                ArrayB1 flux_1g = flux(blitz::Range::all(), ig);

                HDF::Write( node, setname.str(), flux_1g.begin(), flux_1g.end(),
                        dims);
            }

            xs_mesh_->output( node );
            return;
        }


    protected:
        const CoreMesh &mesh_;

        unsigned int n_inner_;
        AngularQuadrature ang_quad_;

        // Boundary condition enumeration
        std::vector<Boundary> bc_type_;

        // Temporary storage for 1-group scalar flux
        ArrayB1 flux_1g_;

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
        template <typename CurrentWorker>
        void sweep_1g( int group ) {
            CurrentWorker cw( coarse_data_, &mesh_ );
            flux_1g_ = 0.0;
            cell_worker_.set_group( group );

            int nx = mesh_.nx();
            int ny = mesh_.ny();
            int nz = mesh_.nz();

            ArrayF x_flux(ny*nz);
            ArrayF y_flux(nx*nz);
            ArrayF z_flux(nx*ny);

            int iang = 0;
            for( auto ang: ang_quad_ ) {
//std::cout << "angle: " << iang << std::endl;
                // Configure the current worker for this angle
                cw.set_octant( iang / ang_quad_.ndir_oct() + 1 );

                cell_worker_.set_angle( iang, ang );

                real_t wgt = ang.weight * HPI;
                real_t ox = ang.ox;
                real_t oy = ang.oy;
                real_t oz = ang.oz;

                // Configure the loop direction. Could template the below for
                // speed at some point.
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
                    cell_worker_.set_z(iz);
                    for( int iy=stty; iy!=stpy; iy+=ydir ) {
                        cell_worker_.set_y(iy);
                        for( int ix=sttx; ix!=stpx; ix+=xdir ) {
                            // Gross. really need an Sn mesh abstraction
                            real_t psi_x = x_flux[ny*iz + iy];
                            real_t psi_y = y_flux[nx*iz + ix];
                            real_t psi_z = z_flux[nx*iy + ix];
//std::cout << psi_x << " " << psi_y << " " << psi_z << std::endl;

                            size_t i = mesh_.coarse_cell( Position( ix, iy, iz ) );

                            //real_t psi = cell_worker_.evaluate( psi_x, psi_y,
                            //        psi_z, q_[i], xstr_[i], i );
                            real_t psi = cell_worker_.evaluate_2d( psi_x, psi_y,
                                    q_[i], xstr_[i], i );

//std::cout << psi << std::endl;
//std::cout << psi_x << " " << psi_y << " " << psi_z << std::endl;
//std::cout << std::endl;

                            x_flux[ny*iz + iy] = psi_x;
                            y_flux[nx*iz + ix] = psi_y;
                            z_flux[nx*iy + ix] = psi_z;

                            flux_1g_(i) += psi*wgt;

                            // Stash currents (or not, depending on the
                            // CurrentWorker template parameter)
                            cw.current_work( x_flux[ny*iz + iy],
                                             y_flux[nx*iz + ix],
                                             z_flux[nx*iy + ix],
                                             i, ang, group );
                        }
                    }
                }

                // store the downwind boundary condition
                bc_out_.set_face(0, iang, Normal::X_NORM, x_flux);
                bc_out_.set_face(0, iang, Normal::Y_NORM, y_flux);
                bc_out_.set_face(0, iang, Normal::Z_NORM, z_flux);
                bc_in_.update( group, iang, bc_out_ );
                iang++;
            }
            // Update the boundary condition
            //bc_in_.update( group, bc_out_ );

            return;
        }

        /**
         * \brief Check the neutron balance in all of the cells of the sweeper
         */
        void check_balance( int group ) const {
            if( !coarse_data_ ) {
                throw EXCEPT("No coarse data. Need it to look at currents.");
            }
            for( size_t icell=0; icell<mesh_.n_pin(); icell++ ) {
                real_t b = 0.0;
    
                // Current
                b -= coarse_data_->current(
                        mesh_.coarse_surf(icell, Surface::EAST), group ) *
                        mesh_.coarse_area( icell, Surface::EAST );
                b -= coarse_data_->current(
                        mesh_.coarse_surf(icell, Surface::NORTH), group ) *
                        mesh_.coarse_area( icell, Surface::NORTH );
                b -= coarse_data_->current(
                        mesh_.coarse_surf(icell, Surface::TOP), group ) *
                        mesh_.coarse_area( icell, Surface::TOP );
                b += coarse_data_->current(
                        mesh_.coarse_surf(icell, Surface::WEST), group ) *
                        mesh_.coarse_area( icell, Surface::WEST );
                b += coarse_data_->current(
                        mesh_.coarse_surf(icell, Surface::SOUTH), group ) *
                        mesh_.coarse_area( icell, Surface::SOUTH );
                b += coarse_data_->current(
                        mesh_.coarse_surf(icell, Surface::BOTTOM), group ) *
                        mesh_.coarse_area( icell, Surface::BOTTOM );
    
                // Source
                b += (*source_)[icell]*vol_[icell];
    
                // Internal removal
                b -= flux_1g_(icell) *
                    (*xs_mesh_)[icell].xsmacrm()[group] * vol_[icell];
    
                std::cout << "Cell balance: " << b << std::endl;
            }
            std::cout << std::endl;
        }
    private:
        Worker cell_worker_;
    };
} }
