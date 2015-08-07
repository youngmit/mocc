#pragma once

#include "xs_mesh.hpp"

namespace mocc {
    class XSMeshHomogenized: public XSMesh {
    public:
        XSMeshHomogenized( const CoreMesh& mesh );
    private:
        const CoreMesh& mesh_;
        XSMeshRegion homogenize_region( int i, int first_reg,
                const Pin& pin ) const;
    };
}
