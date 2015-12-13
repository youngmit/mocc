#pragma once

#include <vector>
#include <memory>

#include <blitz/array.h>

#include "mocc-core/blitz_typedefs.hpp"
#include "mocc-core/eigen_interface.hpp"
#include "mocc-core/global_config.hpp"
#include "mocc-core/mesh.hpp"

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
            flux( (int)mesh.n_pin(), ngroup ),
            old_flux( (int)mesh.n_pin(), ngroup ),
            mesh_(mesh),
            has_data_radial_( false ),
            has_data_axial_( false )
        {
            current = 0.0;
            surface_flux = 0.0;
            flux = 0.0;
            old_flux = 0.0;
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
            current( blitz::Range::all(), group ) = 0.0;
            surface_flux( blitz::Range::all(), group ) = 0.0;
        }

        /**
         * \brief Zero out the data on the radial-normal surfaces for a given
         * group.
         *
         * This is the 2-D version of \ref zero_data(). It zeros out the X- and
         * Y-normal surfaces, but leaves data for the other surfaces untouched.
         */
        void zero_data_radial( int group ) {
            ArrayB1 current_g = current(blitz::Range::all(), group);
            ArrayB1 surface_flux_g = surface_flux(blitz::Range::all(), group);
            for( size_t plane=0; plane<mesh_.nz(); plane++ ) {
                for( auto surf=mesh_.plane_surf_xy_begin(plane);
                        surf!=mesh_.plane_surf_end(plane);
                        ++surf )
                {
                    current_g(surf) = 0.0;
                    surface_flux_g(surf) = 0.0;
                }
            }
            return;
        }

        ArrayB2 current;
        ArrayB2 surface_flux;
        ArrayB2 flux;
        ArrayB2 old_flux;
    private:
        const Mesh &mesh_;
        bool has_data_radial_;
        bool has_data_axial_;
    };

    typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
