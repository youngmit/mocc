#include "transport_sweeper.hpp"

#include <iostream>

namespace mocc {
    real_t TransportSweeper::total_fission( bool old ) const {
        real_t tfis = 0.0;
        const ArrayX& flux = old ? flux_old_: flux_;
        for( auto &xsr: *xs_mesh_ ) {
            for( unsigned int ig=0; ig<ng_; ig++ ) {
                real_t xsnf = xsr.xsmacnf()[ig];
                for (auto &ireg: xsr.reg() ) {
                    tfis += flux(ireg, ig)*vol_(ireg)*xsnf;
                }
            }
        }
        return tfis;
    }

    void TransportSweeper::calc_fission_source( real_t k, 
            ArrayX& fission_source ) const {

        real_t rkeff = 1.0/k;
        fission_source.fill(0.0);
        for( auto &xsr: *xs_mesh_ ) {
            const auto& xsnf = xsr.xsmacnf();
            for( size_t ig=0; ig<xs_mesh_->n_group(); ig++ ) {
                for( auto &ireg: xsr.reg() ) {
                    fission_source(ireg) += rkeff*xsnf[ig]*flux_(ireg, ig);
                }
            }
        }
        return;
    }
}
