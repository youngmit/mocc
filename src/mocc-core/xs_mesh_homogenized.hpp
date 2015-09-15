#pragma once

#include <memory>

#include "eigen_interface.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    class XSMeshHomogenized: public XSMesh {
    public:
        XSMeshHomogenized( const CoreMesh& mesh );
        void update( const ArrayX &flux );
    private:
        const CoreMesh& mesh_;
        XSMeshRegion homogenize_region( int i, int first_reg,
                const Pin& pin ) const;
    };

    typedef std::shared_ptr<XSMeshHomogenized> SP_XSMeshHomogenized_t;
}
