#pragma once

#include "core/source.hpp"

namespace mocc {
    /**
     * This extends the virtual \ref Source type for use as an isotropic source
     * for MoC sweepers.
     */
    class SourceIsotropic : public Source {
    public:
        SourceIsotropic( int nreg, const XSMesh* xs_mesh, const ArrayB2& flux ):
            Source( nreg, xs_mesh, flux),
            q_( nreg )

        {
            return;
        }

        virtual void self_scatter( size_t ig );

        const VectorX& get_transport( int iang ) const {
            return q_;
        }

    protected:
        // The source, including self-scatter. This is stored separately from
        // source_1g_ so that the self_scatter method may be called multiple
        // times without having to completely reconstruct the source. All calls
        // to get_transport() will return a reference to this vector.
        VectorX q_;
    };
}
