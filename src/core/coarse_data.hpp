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

#include <memory>
#include <vector>

#include <blitz/array.h>

#include "core/blitz_typedefs.hpp"
#include "core/eigen_interface.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"

namespace mocc {
    /**
     * CoarseData stores the data needed to do CMFD. Coarse surface currents,
     * fluxes, etc.
     */
    struct CoarseData {
    public:
        CoarseData( const Mesh &mesh, size_t ngroup ):
            current( (int)mesh.n_surf(), ngroup ),
            surface_flux( (int)mesh.n_surf(), ngroup ),
            partial_current( (int)mesh.n_surf(), ngroup ),
            flux( (int)mesh.n_pin(), ngroup ),
            old_flux( (int)mesh.n_pin(), ngroup ),
            n_group_( ngroup ),
            mesh_(mesh),
            has_data_radial_( false ),
            has_data_axial_( false )
        {
            current = 0.0;
            surface_flux = 0.0;
            flux = 0.0;
            old_flux = 0.0;

            //assert( current(blitz::Range::all(), 0).isStorageContiguous() );
            //assert( surface_flux(blitz::Range::all(), 0).
            //        isStorageContiguous() );
            //assert( partial_current(blitz::Range::all(), 0).
            //        isStorageContiguous() );

            return;
        }

        void set_has_radial_data( bool has ) {
            has_data_radial_ = has;
            return;
        }

        void set_has_axial_data( bool has ) {
            has_data_axial_ = has;
            return;
        }

        bool has_axial_data() const {
            return has_data_axial_;
        }

        bool has_radial_data() const {
            return has_data_radial_;
        }

        /**
         * \brief Zero out all of the data associated with the given group.
         *
         * This is typically used immediately before invoking a sweep procedure
         * that will calculate new data.
         *
         * \note This will zero out all of the data, including radial and axial
         * surfaces. It is therefore best suited for use with 3-D sweepers. Most
         * 2-D sweepers will want to use the 2-D version, \ref
         * zero_data_radial().
         */
        void zero_data( int group ) {
            assert(group < n_group_ );
            // need this to disambiguate the operator= for partial currents.
            const std::array<real_t, 2> zero = {0.0, 0.0};

            current( blitz::Range::all(), group ) = 0.0;
            surface_flux( blitz::Range::all(), group ) = 0.0;
            partial_current( blitz::Range::all(), group ) = zero;
        }

        /**
         * \brief Zero out the data on the radial-normal surfaces for a given
         * group.
         *
         * This is the 2-D version of \ref zero_data(). It zeros out the X- and
         * Y-normal surfaces, but leaves data for the other surfaces untouched.
         */
        void zero_data_radial( int group ) {
            assert(group < n_group_ );

            // need this to disambiguate the operator= for partial currents.
            const std::array<real_t, 2> zero = {0.0, 0.0};

            auto current_g = current(blitz::Range::all(), group);
            auto surface_flux_g = surface_flux(blitz::Range::all(), group);
            auto partial_g = partial_current(blitz::Range::all(), group);
            for( size_t plane=0; plane<mesh_.nz(); plane++ ) {
                for( auto surf=mesh_.plane_surf_xy_begin(plane);
                        surf!=mesh_.plane_surf_end(plane);
                        ++surf )
                {
                    current_g(surf) = 0.0;
                    surface_flux_g(surf) = 0.0;
                    partial_g(surf) = zero;
                }
            }
            return;
        }

        ArrayB2 current;
        ArrayB2 surface_flux;
        blitz::Array<std::array<real_t, 2>, 2> partial_current;
        ArrayB2 flux;
        ArrayB2 old_flux;
    private:
        int n_group_;
        const Mesh &mesh_;
        bool has_data_radial_;
        bool has_data_axial_;
    };

    typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
