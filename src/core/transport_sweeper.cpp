/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

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
            for( int ig=0; ig<n_group_; ig++ ) {
                real_t xsnf = xsr.xsmacnf(ig);
                for (auto &ireg: xsr.reg() ) {
                    tfis += flux((int)ireg, (int)ig)*vol_[ireg]*xsnf;
                }
            }
        }
        return tfis;
    }

    void TransportSweeper::calc_fission_source( real_t k,
            ArrayB1& fission_source ) const {

        real_t rkeff = 1.0/k;
        fission_source = 0.0;
        for( auto &xsr: *xs_mesh_ ) {
            const auto &xsnf = xsr.xsmacnf();
            for( int ig=0; ig<(int)n_group_; ig++ ) {
                for( auto &ireg: xsr.reg() ) {
                    fission_source(ireg) += rkeff * xsnf[ig] *
                        flux_(ireg, ig);
                }
            }
        }

        return;
    }

    ArrayB3 TransportSweeper::pin_powers() const {
        assert(n_reg_ == (int)core_mesh_->n_reg());
        ArrayB3 powers(core_mesh_->nz(), core_mesh_->ny(), core_mesh_->nx());
        powers = 0.0;

        // This isnt the most efficient way to do this, memory-wise, but its
        // quick and simple. Calculate volume x flux x kappa-fission for all
        // flat source regions, then reduce to the pin mesh.
        ArrayB1 fsr_pow(n_reg_);
        fsr_pow = 0.0;
        for( const auto &xsr: *xs_mesh_ ) {
            for( int ig=0; ig<n_group_; ig++ ) {
                for( auto ireg: xsr.reg() ) {
                    fsr_pow(ireg) += flux_(ireg, ig) * xsr.xsmackf(ig) *
                        vol_[ireg];
                }
            }
        }


        int ipin = 0;
        int ireg = 0;
        real_t tot_pow = 0.0;
        for( const auto pin: *core_mesh_ ) {
            if( pin->is_fuel() ) {
                const PinMesh &pm = pin->mesh();
                Position pos = core_mesh_->pin_position(ipin);
                for( int ir=0; ir<pm.n_reg(); ir++ ) {
                    tot_pow += fsr_pow(ireg);
                    powers(pos.z, pos.y, pos.x) += fsr_pow(ireg);
                    ireg++;
                }
            }
            ipin++;
        }

        // Normalize!
        tot_pow = powers.size()*core_mesh_->n_fuel_2d()/tot_pow;
        for( auto &v: powers ) {
            v *= tot_pow;
        }


        return powers;
    }

    ArrayB2 TransportSweeper::get_pin_flux() const {
        assert( core_mesh_ );
        ArrayB2 flux( core_mesh_->n_pin(), n_group_ );

        auto flux_it = flux.begin();

        for( int ig=0; ig<n_group_; ig++ ) {
            ArrayB1 flux_1g( flux(blitz::Range::all(), ig) );
            this->get_pin_flux_1g( ig, flux_1g );
        }

        return flux;
    }

    real_t TransportSweeper::flux_residual() const {
        real_t r = 0.0;
        auto it = flux_.begin();
        auto it_old = flux_old_.begin();
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
