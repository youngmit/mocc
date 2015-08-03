#include "transport_sweeper.hpp"

namespace mocc {
    float_t TransportSweeper::total_fission( bool old ) const {
        float_t tfis = 0.0;
        const ArrayX& flux = old ? flux_old_: flux_;
        for( auto &xsr: xs_mesh_ ) {
            for( int ig=0; ig<ng_; ig++ ) {
                float_t xsnf = xsr.xsmacnf()[ig];
                for (auto &ireg: xsr.reg() ) {
                    tfis += flux(ireg, ig)*vol_(ireg)*xsnf;
                }
            }
        }
        return tfis;
    }
}
