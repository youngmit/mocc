#include "sn_sweeper_dd.hpp"

namespace mocc {
    void SnSweeper_DD::sweep( int group ) {
        // Store the transport cross section somewhere useful
        for( auto &xsr: *xs_mesh_ ) {
            real_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_[ireg] = xstr;
            }
        }
        
        flux_1g_ = flux_[ std::slice(n_reg_*group, n_reg_, 1) ];

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            // Set the source (add upscatter and divide by 4PI)
            source_->self_scatter( group, flux_1g_, q_ );

            if( inner == n_inner_-1 && coarse_data_ ) {
                // Wipe out the existing currents
                coarse_data_->current.col( group ) = 0.0;
                this->sweep_1g<sn::Current, CellWorker_DD>( group, 
                        cell_worker_ );
            } else {
                this->sweep_1g<sn::NoCurrent, CellWorker_DD>( group, 
                        cell_worker_ );
            }
        }
        flux_[ std::slice(n_reg_*group, n_reg_, 1) ] = flux_1g_;

        return;
    }
}