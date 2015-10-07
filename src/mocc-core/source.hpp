#pragma once

#include <iostream>

#include "eigen_interface.hpp"

#include "xs_mesh.hpp"

namespace mocc {
    class Source {
    public:
        Source( int nreg, const XSMesh* xs_mesh, const ArrayX& flux );

        virtual ~Source(){
        }

        // Initialize the source with the external source, if it exists, then
        // add the groups contribution from the multi-group fission source
        virtual void fission( const ArrayX& fs, int ig );

        // Add the contribution from in-scattering from other groups. At some
        // point, ill play with upscattering iterations, but for now KISS.
        virtual void in_scatter( size_t ig );

        // Add a contribution due to self-scatter within the current group,
        // returning the final source. This is usually called several times by a
        // sweeper in its "inner" iterations, and therefore does not mutate the
        // interal representation of the source, but instead returns the result
        // to the caller through the qbar argument.
        virtual void self_scatter( size_t ig, ArrayX& flux_1g, 
                ArrayX& qbar ) const;

        // Return a pointer to the source
        const real_t* get() const {
            return source_1g_.data();
        }
    protected:
        const XSMesh *xs_mesh_;
        size_t ng_;
        // This is true if an external source has been specified. For now it's
        // initialized false.
        bool has_external_;

        // Single-group source
        ArrayX source_1g_;

        // Reference to the MG flux variable. Need this to do scattering
        // contributions, etc.
        const ArrayX& flux_;
    };

    typedef std::shared_ptr<Source> SP_Source_t;
    typedef std::unique_ptr<Source> UP_Source_t;
}
