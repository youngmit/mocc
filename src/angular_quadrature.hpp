#pragma once

#include <vector>

#include "global_config.hpp"
#include "angle.hpp"

namespace mocc {

    enum quad_t {
        SN // Level-symmetric
    };

    class AngularQuadrature {
    public:
        AngularQuadrature( quad_t type, int azi_order, int pol_order=0 );
    private:
        const quad_t type_;
        std::vector<Angle> angles_;
    };
}
