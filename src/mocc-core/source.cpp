#include "source.hpp"

#include <algorithm>

#include "mocc-core/error.hpp"
#include "mocc-core/constants.hpp"

using std::endl;
using std::cout; 

namespace mocc {
    Source::Source( int nreg, const XSMesh *xs_mesh, const ArrayF& flux ):
        xs_mesh_( xs_mesh ),
        n_group_( xs_mesh->n_group() ),
        has_external_( false ),
        source_1g_( nreg ),
        flux_( flux ),
        n_reg_( flux.size()/n_group_ )
    {
        assert(n_reg_*n_group_ == flux_.size() );
        return;
    }

    // Multiply the group-independent fission source by \c chi[ig] to get the
    // fission source into the current group. If an external source is defines,
    // start with that.
    void Source::fission( const ArrayF& fs, int ig ) {
        assert(fs.size() == n_reg_);
        if( has_external_ ) {
            throw EXCEPT( "No support for external sources yet." );
        } else {
            source_1g_ = 0.0;
        }

        for( auto &xsr: *xs_mesh_ ) {
            real_t xsch = xsr.xsmacch()[ig];
            for( auto &ireg: xsr.reg() ) {
                source_1g_[ireg] +=  xsch * fs[ireg];
            }
        }
        return;
    }

    // Compute the contribution to the source from inscattering from other
    // groups
    void Source::in_scatter( size_t ig ) {
        for( auto &xsr: *xs_mesh_ ) {
            const ScatteringRow& scat_row = xsr.xsmacsc().to(ig);
            size_t min_g = scat_row.min_g;
            size_t igg = min_g;
            for( auto sc: scat_row ) {
                // Dont add a contribution for self-scatter. TODO: it might be a
                // good idea to remove self-scatter from the scattering matrix
                // and store it separately. May also benefit the self-scatter
                // routine to have less indirection.
                if( igg != ig ) {
                    for( auto &ireg: xsr.reg() ) {
                        real_t scat_src = sc*flux_[ireg + n_reg_*igg];
                        source_1g_[ireg] += scat_src;
                    }
                }
                igg++;
            }
        }
        return;
    }

    // This can get away with being const, since we are actually returning the
    // source to the caller. Nothing should get touched internally
    void Source::self_scatter( size_t ig, ArrayF& flux_1g, 
            ArrayF& qbar ) const {
        for( auto &xsr: *xs_mesh_ ) {
            const ScatteringRow& scat_row = xsr.xsmacsc().to(ig);
            real_t xssc = scat_row.from[ig-scat_row.min_g];
            real_t r_fpi_tr = 1.0/(xsr.xsmactr()[ig]*FPI);
            for ( auto &ireg: xsr.reg() ) {
                qbar[ireg] = ( source_1g_[ireg] + flux_1g[ireg]*xssc ) * 
                    r_fpi_tr;
            }
        }

        // Check to make sure that the source is positive
        /*if( std::any_of(qbar.begin(), qbar.end(), [](real_t v){ return v < 0.0; }) ) {
            throw EXCEPT("Negative source!");
        }*/

        return;
    }
}
