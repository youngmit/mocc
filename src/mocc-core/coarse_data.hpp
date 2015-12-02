#pragma once

#include <vector>
#include <memory>

#include <blitz/array.h>

#include "arrays.hpp"
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
            has_data_( false )
        {
            current = 0.0;
            flux = 0.0;
            old_flux = 0.0;
            return;
        }

        void set_has_data( bool has ) {
            has_data_ = has;
            return;
        }

        bool has_data() const {
            return has_data_;
        }

        ArrayB2 current;
        ArrayB2 flux;
        ArrayB2 old_flux;
    private:
        bool has_data_;
    };

    typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
