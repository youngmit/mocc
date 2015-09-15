#pragma once

#include <vector>
#include <memory>

#include "arrays.hpp"
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
            old_flux( nreg, ngroup )
        {
            return;
        }

        ArrayX current;
        ArrayX flux;
        ArrayX old_flux;
    };

    typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
