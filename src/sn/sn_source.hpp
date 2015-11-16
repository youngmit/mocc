#pragma once

#include <iostream>

#include "source.hpp"

namespace mocc {
    /**
     * Slight modification to the original \ref Source type to avoid scaling the
     * source by the transport cross section, which is unneccesary for the Sn
     * sweepers. Could do some template magic at some point instead.
     */
    class SnSource: public Source {
    public:
        SnSource( int nreg, const XSMesh *xs_mesh, const ArrayF& flux ):
            Source( nreg, xs_mesh, flux )
        {
            return;
        }

        ~SnSource() {
        }

        /**
         * Replaces the standard \ref Source::self_scatter() method with one that
         * does not divide the source by the transport cross section, which is
         * only needed for the MoC sweeper.
         */
        void self_scatter( size_t ig, ArrayF& flux_1g,
                ArrayF& qbar ) const
        {
            for( auto &xsr: *xs_mesh_ ) {
                const ScatteringRow& scat_row = xsr.xsmacsc().to(ig);
                real_t xssc = scat_row.from[ig-scat_row.min_g];
                for ( auto &ireg: xsr.reg() ) {
                    qbar[ireg] = ( source_1g_[ireg] + flux_1g[ireg]*xssc ) *
                        RFPI;
                }
            }
            return;
        }
    };
}
