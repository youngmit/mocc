#include "angular_quadrature.hpp"

#include <vector>

#include "error.hpp"
#include "level_symmetric.hpp"


namespace mocc {
    AngularQuadrature::AngularQuadrature( quad_t type, int azi_order, 
                                          int pol_order ):
        type_(type){
        
        switch(type) {
            case SN:
                angles_ = GenSn( azi_order );
                break;
        }

        return;
    }
}
