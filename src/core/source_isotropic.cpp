#include "core/source.hpp"
#include "core/source_isotropic.hpp"

namespace mocc {
    void SourceIsotropic::self_scatter( size_t ig )
    {
        // Take a slice reference for this group's flux
        const ArrayB1 flux_1g = flux_(blitz::Range::all(), ig);
        for( auto &xsr: *xs_mesh_ ) {
            const ScatteringRow& scat_row = xsr.xsmacsc().to(ig);
            real_t xssc = scat_row.from[ig-scat_row.min_g];
            real_t r_fpi_tr = 1.0/(xsr.xsmactr()[ig]*FPI);
            for ( auto &ireg: xsr.reg() ) {
                q_[ireg] = ( source_1g_[ireg] + flux_1g((int)ireg)*xssc ) *
                    r_fpi_tr;
            }
        }

        // Check to make sure that the source is positive
        bool any = false;
        for( int i=0; i<q_.size(); i++ ) {
            if(q_[i] < 0.0 ) {
                any = true;
            }
        }
        if( any ) {
          //  throw EXCEPT("Negative source!");
        }

        return;

    }
}
