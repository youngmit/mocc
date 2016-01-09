#pragma once

#include "mocc-core/source.hpp"

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

        virtual void self_scatter( size_t ig, const ArrayB1 &flux_1g );

        const VectorX& get_transport( int iang ) const {
            return q_;
        }

    private:
        VectorX q_;
    };
}
