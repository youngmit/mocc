#pragma once

#include "source.hpp"

namespace mocc {
    // Slight modification to the original Source type to avoid scaling the
    // source by the transport cross section, which is unneccesary for the Sn
    // sweepers. Could do some template magic at some point, instead.
    class SnSource: public Source {
    public:
        SnSource( int nreg, const XSMesh *xs_mesh, const ArrayX& flux ):
            Source( nreg, xs_mesh, flux )
        {
            return;
        }

        void self_scatter( unsigned int ig, ArrayX& flux_1g, 
                ArrayX& qbar ) const 
        {
            for( auto &xsr: *xs_mesh_ ) {
                const ScatRow& scat_row = xsr.xsmacsc().to(ig);
                float_t xssc = scat_row.from[ig-scat_row.min_g];
                for ( auto &ireg: xsr.reg() ) {
                    qbar(ireg) = ( source_1g_(ireg) + flux_1g(ireg)*xssc ) * 
                        RFPI;
                    //qbar(ireg) = ( source_1g_(ireg) ) * RFPI;
                }
            }

            return;
        }
    };
}
