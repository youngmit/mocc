#pragma once

#include "mocc-core/angular_quadrature.hpp"
#include "mocc-core/transport_sweeper.hpp"
#include "mocc-core/utils.hpp"

#include "sn/sn_boundary.hpp"
#include "sn/cell_worker.hpp"
#include "sn/sn_source.hpp"

namespace mocc { namespace sn {
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node &input, const CoreMesh& mesh ):
            TransportSweeper( input ),
            mesh_( mesh ),
            bc_type_( mesh.boundary() ),
            flux_1g_( mesh_.n_pin() ),
            xstr_( mesh.n_pin() ),
            q_( mesh_.n_pin() ),
            bc_in_( mesh.mat_lib().n_group(), ang_quad_, mesh_ ),
            bc_out_( 1, ang_quad_, mesh_ )
        {
            return;
        }

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

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return std::static_pointer_cast<XSMeshHomogenized>( xs_mesh_ );
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

    };

    typedef std::shared_ptr<SnSweeper> SP_SnSweeper_t;
    typedef std::unique_ptr<SnSweeper> UP_SnSweeper_t;
}}