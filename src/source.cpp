#include "source.hpp"

#include "error.hpp"

namespace mocc {
    Source::Source() {
        return;
    }

    Source::Source( int nreg, const XSMesh& xs_mesh ):
        xs_mesh_(&xs_mesh),
        ng_(xs_mesh.n_grp()),
        has_external_(false),
        source_1g_(nreg, 1)
    {
        return;
    }

    // Multiply the group-independent fission source by chi[ig] to get the
    // fission source into the current group. If an external source is defines,
    // start with that.
    void Source::fission( const MatrixX& fs, int ig ) {       
        if( has_external_ ) {
            Error( "No support for external sources yet." );
        } else {
            source_1g_.fill(0.0);
        }

    }
}
