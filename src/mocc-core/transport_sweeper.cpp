#include "transport_sweeper.hpp"

#include <cmath>
#include <iostream>

using std::endl;
using std::cout;

namespace mocc {
    real_t TransportSweeper::total_fission( bool old ) const {
        real_t tfis = 0.0;
        const auto &flux = old ? flux_old_: flux_;
        for( auto &xsr: *xs_mesh_ ) {
            for( unsigned int ig=0; ig<n_group_; ig++ ) {
                real_t xsnf = xsr.xsmacnf()[ig];
                for (auto &ireg: xsr.reg() ) {
                    tfis += flux((int)ireg, (int)ig)*vol_[ireg]*xsnf;
                }
            }
        }
        return tfis;
    }

    void TransportSweeper::calc_fission_source( real_t k,
            ArrayF& fission_source ) const {

        real_t rkeff = 1.0/k;
        fission_source = 0.0;
        for( auto &xsr: *xs_mesh_ ) {
            const auto &xsnf = xsr.xsmacnf();
            for( size_t ig=0; ig<n_group_; ig++ ) {
                for( auto &ireg: xsr.reg() ) {
                    fission_source[ireg] += rkeff * xsnf[ig] *
                        flux_((int)ireg, (int)ig);
                }
            }
        }
        return;
    }

    ArrayB2 TransportSweeper::get_pin_flux() const {
        assert( core_mesh_ );
        ArrayB2 flux( core_mesh_->n_pin(), n_group_ );

        auto flux_it = flux.begin();

        for( size_t ig=0; ig<n_group_; ig++ ) {
            ArrayB1 flux_1g( flux(blitz::Range::all(), ig) );
            this->get_pin_flux_1g( ig, flux_1g );
        }

        return flux;
    }

    real_t TransportSweeper::flux_residual() const {
        real_t r = 0.0;
        auto it = flux_.begin();
        auto it_old = flux_.begin();
        auto end = flux_.end();
        while( it != end ) {
            real_t e = *it - *it_old;
            r += e*e;
            ++it;
            ++it_old;
        }
        return std::sqrt(r);
    }


}
