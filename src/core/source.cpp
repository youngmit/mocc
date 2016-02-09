#include "source.hpp"

#include <algorithm>

#include "core/error.hpp"
#include "core/constants.hpp"

namespace mocc {
    Source::Source( int nreg, const XSMesh *xs_mesh, const ArrayB2& flux ):
        xs_mesh_( xs_mesh ),
        n_group_( xs_mesh->n_group() ),
        has_external_( false ),
        scale_transport_( false ),
        source_1g_( nreg ),
        flux_( flux ),
        n_reg_( flux.size()/n_group_ )
    {
        assert( n_reg_*n_group_ == flux_.size() );
        state_.reset();
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

        state_.reset();
        return;
    }

    // Multiply the group-independent fission source by \c chi[ig] to get the
    // fission source into the current group. If an external source is defines,
    // start with that.
    void Source::fission( const ArrayB1& fs, int ig ) {
        assert(fs.size() == n_reg_);
        assert(!state_.has_fission);
        assert(!state_.is_scaled);

        for( auto &xsr: *xs_mesh_ ) {
            real_t xsch = xsr.xsmacch(ig);
            for( const int &ireg: xsr.reg() ) {
                source_1g_[ireg] +=  xsch * fs(ireg);
            }
        }

        state_.has_fission = true;

        return;
    }

    /**
     * \brief Compute the contribution to the source from inscattering from
     * other groups.
     */
    void Source::in_scatter( size_t ig ) {
        assert(!state_.has_inscatter);
        assert(!state_.is_scaled);
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


    void Source::auxiliary( const ArrayB1 &aux ) {
        assert( source_1g_.size() == (int)aux.size() );
        assert( !state_.is_scaled );
        for( int i=0; i<(int)n_reg_; i++) {
            source_1g_[i] += aux(i);
        }
    }

    void Source::add_external( const pugi::xml_node &input ) {
        if( input.attribute("file").empty() ) {
            // Nothing to do here
            return;
        }

        std::string srcfname = input.attribute("file").value();
        HDF::H5File srcfile( srcfname, "r" );
        VecF src;
        VecI dims;
        HDF::Read( srcfile.get(), "/source", src, dims );

        if( dims[0] != (int)n_group_ ) {
            throw EXCEPT("Wrong group dimensions for source");
        }
        if( dims[1] != (int)n_reg_ ) {
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
