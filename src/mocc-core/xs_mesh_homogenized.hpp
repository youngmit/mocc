#pragma once

#include <memory>

#include "eigen_interface.hpp"
#include "h5file.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    class XSMeshHomogenized: public XSMesh {
    public:
        XSMeshHomogenized( const CoreMesh& mesh );

        /**
         * Update homogenized cross sections using the passed flux array
         */
        void update( const ArrayX &flux );

        /**
         * Generate output of important cross sections on the homogenized mesh
         */
        void output( H5File &file ) const;
    private:
        const CoreMesh& mesh_;
        
        /** 
        * \brief Return an XSMeshRegion containing homogenized cross sections
        * from a pin cell. No flux wieghting is performed, only volume
        * weighting.
        */
        XSMeshRegion homogenize_region( int i, const Pin& pin ) const;
        
        /** 
        * \brief Return an XSMeshRegion containing homogenized cross sections
        * from a pin cell. Use the passed scalar flux to perform flux-volume
        * weighting.
        */
        XSMeshRegion homogenize_region( int i, int first_reg, const Pin& pin,
                const ArrayX &flux ) const;
    };

    typedef std::shared_ptr<XSMeshHomogenized> SP_XSMeshHomogenized_t;
}
