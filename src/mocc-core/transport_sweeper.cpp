#include "transport_sweeper.hpp"

#include <iostream>

namespace mocc {
    float_t TransportSweeper::total_fission( bool old ) const {
        float_t tfis = 0.0;
        const ArrayX& flux = old ? flux_old_: flux_;
        for( auto &xsr: *xs_mesh_ ) {
            for( int ig=0; ig<ng_; ig++ ) {
                float_t xsnf = xsr.xsmacnf()[ig];
                for (auto &ireg: xsr.reg() ) {
                    tfis += flux(ireg, ig)*vol_(ireg)*xsnf;
                }
            }
        }
        return tfis;
    }

    void TransportSweeper::calc_fission_source( float_t k, 
            ArrayX& fission_source ) const {

        float_t rkeff = 1.0/k;
        fission_source.fill(0.0);
        for( auto &xsr: *xs_mesh_ ) {
            const auto& xsnf = xsr.xsmacnf();
            for( unsigned int ig=0; ig<xs_mesh_->n_grp(); ig++ ) {
                for( auto &ireg: xsr.reg() ) {
                    fission_source(ireg) += rkeff*xsnf[ig]*flux_(ireg, ig);
                }
            }
        }
        return;
    }
}
