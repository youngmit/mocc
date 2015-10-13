#include "transport_sweeper.hpp"

#include <iostream>

namespace mocc {
    real_t TransportSweeper::total_fission( bool old ) const {
        real_t tfis = 0.0;
        const ArrayX& flux = old ? flux_old_: flux_;
        for( auto &xsr: *xs_mesh_ ) {
            for( unsigned int ig=0; ig<n_group_; ig++ ) {
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

    VecF TransportSweeper::get_pin_flux() const {
        VecF flux( core_mesh_->n_pin()*n_group_, 0.0 );

        auto flux_it = flux.begin();

        VecF flux_1g;
        for( int ig=0; ig<n_group_; ig++ ) {
            this->get_pin_flux_1g( ig, flux_1g );
            flux_it = std::copy( flux_1g.begin(), flux_1g.end(), flux_it );
        }

        return flux;
    }
}
