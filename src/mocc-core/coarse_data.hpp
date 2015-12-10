#pragma once

#include <vector>
#include <memory>

#include <blitz/array.h>

#include "blitz_typedefs.hpp"
#include "eigen_interface.hpp"
#include "global_config.hpp"

namespace mocc {
    /**
     * CoarseData stores the data needed to do CMFD. Coarse surface currents,
     * fluxes, etc.
     */
    struct CoarseData {
    public:
        CoarseData( size_t nreg, size_t nsurf, size_t ngroup ):
            current( nsurf, ngroup ),
            flux( nreg, ngroup ),
            old_flux( nreg, ngroup ),
            has_data_radial_( false ),
            has_data_axial_( false )
        {
            current = 0.0;
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

        ArrayB2 current;
        ArrayB2 flux;
        ArrayB2 old_flux;
    private:
        bool has_data_radial_;
        bool has_data_axial_;
    };

    typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
