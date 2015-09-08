#pragma once

#include <vector>
#include <memory>

#include "global_config.hpp"

namespace mocc {
    /**
     * CoarseData stores the data needed to do CMFD. Coarse surface currents,
     * fluxes, etc.
     */
    class CoarseData {
    public:
        std::vector<VecF> get_current() {
            return current_;
        }
        
        VecF get_gurrent( int group ) {
            return current_[group];
        }

        std::vector<VecF> get_flux() {
            return flux_;
        }

        VecF get_flux( int group ) {
            return flux_[group];
        }

    private:
        std::vector<VecF> current_;
        std::vector<VecF> flux_;
        std::vector<VecF> old_flux_;
    };

    typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
