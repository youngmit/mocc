#include "source.hpp"

#include "error.hpp"
#include "constants.hpp"

namespace mocc {
    Source::Source( int nreg, const XSMesh *xs_mesh, const ArrayX& flux ):
        xs_mesh_(xs_mesh),
        ng_(xs_mesh->n_grp()),
        has_external_(false),
        source_1g_(nreg, 1),
        flux_(flux)
    {
        return;
    }

    // Multiply the group-independent fission source by chi[ig] to get the
    // fission source into the current group. If an external source is defines,
    // start with that.
    void Source::fission( const ArrayX& fs, int ig ) {
        if( has_external_ ) {
            throw EXCEPT( "No support for external sources yet." );
        } else {
            source_1g_.fill(0.0);
        }

        for( auto &xsr: *xs_mesh_ ) {
            float_t xsch = xsr.xsmacch()[ig];
            for( auto &ireg: xsr.reg() ) {
                source_1g_(ireg) +=  xsch*fs(ireg);
            }
        }
        return;
    }

    // Compute the contribution to the source from inscattering from other
    // groups
    void Source::in_scatter( unsigned int ig ) {
        for( auto &xsr: *xs_mesh_ ) {
            const ScatRow& scat_row = xsr.xsmacsc().to(ig);
            unsigned int min_g = scat_row.min_g;
            int igg = min_g;
            for( auto sc: scat_row ) {
                // Dont add a contribution for self-scatter. TODO: it might be a
                // good idea to remove self-scatter from the scattering matrix
                // and store it separately. May also benefit the self-scatter
                // routine to have less indirection.
                if( igg != ig ) {
                    for( auto &ireg: xsr.reg() ) {
                        float_t scat_src = sc*flux_(ireg, igg);
                        source_1g_(ireg) += scat_src;
                    }
                }
                igg++;
            }
        }
        return;
    }

    // This can get away with being const, since we are actually returning the
    // source to the caller. Nothing should get touched internally
    void Source::self_scatter( unsigned int ig, ArrayX& flux_1g, 
            ArrayX& qbar ) const {
        for( auto &xsr: *xs_mesh_ ) {
            const ScatRow& scat_row = xsr.xsmacsc().to(ig);
            float_t xssc = scat_row.from[ig-scat_row.min_g];
            float_t r_fpi_tr = 1.0/(xsr.xsmactr()[ig]*FPI);
            for ( auto &ireg: xsr.reg() ) {
                qbar(ireg) = ( source_1g_(ireg) + flux_1g(ireg)*xssc ) * 
                    r_fpi_tr;
            }
        }

        return;
    }
}
