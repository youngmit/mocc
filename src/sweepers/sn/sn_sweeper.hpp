/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include "core/angular_quadrature.hpp"
#include "core/boundary_condition.hpp"
#include "core/pugifwd.hpp"
#include "core/timers.hpp"
#include "core/transport_sweeper.hpp"
#include "core/utils.hpp"

#include "cell_worker.hpp"

namespace mocc { namespace sn {
    class SnSweeper: public TransportSweeper {
    public:
        SnSweeper( const pugi::xml_node &input, const CoreMesh &mesh );

        void initialize() {
            flux_ = 1.0;
            flux_old_ = 1.0;
            bc_in_.initialize_scalar(1.0/FPI);

            return;
        }

        void get_pin_flux_1g( int ig, ArrayB1& flux ) const {
            assert( (int)flux.size() == n_reg_ );

            flux = flux_(blitz::Range::all(), ig);

            return;
        }

        ArrayB3 pin_powers() const;

        /**
         * Just copy the flux across, since no homogenization is necessary.
         */
        real_t set_pin_flux_1g( int group, const ArrayB1 &pin_flux ) {
            assert( (int)pin_flux.size() == n_reg_ );

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

        virtual void output( H5Node &node ) const;

    protected:
        Timer &timer_;
        Timer &timer_init_;
        Timer &timer_sweep_;
        const CoreMesh &mesh_;

        unsigned int n_inner_;

        // Boundary condition enumeration
        std::array<Boundary, 6> bc_type_;

        // One-group slice of flux_. Should be default-constructed, and assigned
        // slices using .reference()
        ArrayB1 flux_1g_;

        // Temporary storage of the current-group transport cross section
        ArrayF xstr_;

        // Incomming boundary condition
        BoundaryCondition bc_in_;

        // Outgoing boundary condition. Only difined for one group
        BoundaryCondition bc_out_;

        // Gauss-Seidel BC update?
        bool gs_boundary_;
        
        // Protected methods
        /** 
         * \brief Grab data (XS, etc.) from one or more external files
         */
        void add_data( const pugi::xml_node &input );

        /**
         * \brief Check the neutron balance in all of the cells of the sweeper
         */
        void check_balance( int group ) const;

    private:

    };

    typedef std::shared_ptr<SnSweeper> SP_SnSweeper_t;
    typedef std::unique_ptr<SnSweeper> UP_SnSweeper_t;
}}
