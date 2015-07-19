#pragma once

#include "eigen_interface.hpp"

#include "xs_mesh.hpp"

namespace mocc {
    class Source {
    public:
        Source( int nreg, const XSMesh& xs_mesh, const MatrixX& flux );

        // Initialize the source with the external source, if it exists, then
        // add the groups contribution from the multi-group fission source
        void fission( const MatrixX& fs, int ig );

        // Add the contribution from in-scattering from other groups. At some
        // point, ill play with upscattering iterations, but for now KISS.
        void in_scatter( unsigned int ig );

        // Add a contribution due to self-scatter within the current group,
        // returning the final source. This is usually called several times by a
        // sweeper in its "inner" iterations, and therefore does not mutate the
        // interal representation of the source, but instead returns the result
        // to the caller through the qbar argument.
        void self_scatter( unsigned int ig, MatrixX& qbar ) const;

        // Return a pointer to the source
        const float_t* get() const {
            return source_1g_.data();
        }
    private:
        const XSMesh& xs_mesh_;
        unsigned int ng_;
        // This is true if an external source has been specified. For now it's
        // initialized false.
        bool has_external_;

        // Single-group source
        MatrixX source_1g_;

        // Reference to the MG flux variable. Need this to do scattering
        // contributions, etc.
        const MatrixX& flux_;
    };

    typedef std::shared_ptr<Source> SP_Source_t;
}
