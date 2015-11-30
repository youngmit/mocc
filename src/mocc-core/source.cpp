#include "source.hpp"

#include <algorithm>

#include "mocc-core/error.hpp"
#include "mocc-core/constants.hpp"

using std::endl;
using std::cout;

namespace mocc {
    Source::Source( int nreg, const XSMesh *xs_mesh, const ArrayB2& flux ):
        xs_mesh_( xs_mesh ),
        n_group_( xs_mesh->n_group() ),
        has_external_( false ),
        source_1g_( nreg ),
        flux_( flux ),
        n_reg_( flux.size()/n_group_ )
    {
        assert( n_reg_*n_group_ == flux_.size() );
        return;
    }

    void Source::initialize_group( int ig ) {
        if( has_external_ ) {
            for( int ireg=0; ireg<(int)n_reg_; ireg++ ) {
                source_1g_[ireg] = external_source_(ireg, ig);
            }
        } else {
            source_1g_.fill(0.0);
        }
        return;
    }

    // Multiply the group-independent fission source by \c chi[ig] to get the
    // fission source into the current group. If an external source is defines,
    // start with that.
    void Source::fission( const ArrayF& fs, int ig ) {
        assert(fs.size() == n_reg_);

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
            if( xsr.reg().size() == 0) {
                continue;
            }
            const ScatteringRow& scat_row = xsr.xsmacsc().to(ig);
            size_t min_g = scat_row.min_g;
            int igg = min_g;
            for( auto sc: scat_row ) {
                // Dont add a contribution for self-scatter. TODO: it might be a
                // good idea to remove self-scatter from the scattering matrix
                // and store it separately. May also benefit the self-scatter
                // routine to have less indirection.
                if( igg != (int)ig ) {
                    for( auto &ireg: xsr.reg() ) {
                        real_t scat_src = sc*flux_((int)ireg, igg);
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
    void Source::self_scatter( size_t ig, const ArrayB1& flux_1g,
            ArrayF& qbar ) const {
        for( auto &xsr: *xs_mesh_ ) {
            const ScatteringRow& scat_row = xsr.xsmacsc().to(ig);
            real_t xssc = scat_row.from[ig-scat_row.min_g];
            real_t r_fpi_tr = 1.0/(xsr.xsmactr()[ig]*FPI);
            for ( auto &ireg: xsr.reg() ) {
                qbar[ireg] = ( source_1g_[ireg] + flux_1g((int)ireg)*xssc ) *
                    r_fpi_tr;
            }
        }

        // Check to make sure that the source is positive
        bool any = false;
        for( size_t i=0; i<qbar.size(); i++ ) {
            if(qbar[i] < 0.0 ) {
                any = true;
            }
        }
        if( any ) {
          //  throw EXCEPT("Negative source!");
        }

        return;
    }

    void Source::auxiliary( const ArrayB1 &aux ) {
        assert( source_1g_.size() == (int)aux.size() );
        for( int i=0; i<(int)n_reg_; i++) {
            source_1g_[i] += aux(i);
        }
    }

    void Source::add_external( const pugi::xml_node &input ) {
        // Actual source specification. Pretty limited for now.
        if( input.empty() ) {
            throw EXCEPT("Standalone FSS must supply a <source> "
                    "specification");
        }
        std::string srcfname = input.attribute("file").value();
        HDF::H5File srcfile( srcfname, "r" );
        VecF src;
        VecI dims;
        HDF::Read( srcfile.get(), "/source", src, dims );

        if( dims[0] != n_group_ ) {
            throw EXCEPT("Wrong group dimensions for source");
        }
        if( dims[1] != n_reg_ ) {
            throw EXCEPT("Wrong regions dimensions for source");
        }

        external_source_.resize( n_reg_, n_group_ );
        // Do a 1:1 copy. We are assuming that the data in the vector is laid
        // out in the same order as the ArrayB2.
        auto in_it = src.begin();
        auto it = external_source_.begin();
        auto end = src.end();
        while( in_it != end ) {
            *it = *in_it;
            ++it;
            ++in_it;
        }

        has_external_ = true;
    }
}
