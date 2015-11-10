#pragma once

#include <iostream>

#include "eigen_interface.hpp"

#include "xs_mesh.hpp"

namespace mocc {
    class Source {
    public:
        Source( int nreg, const XSMesh* xs_mesh, const ArrayF& flux );

        virtual ~Source(){
        }

        /**
         * Initialize the source with the external source, if it exists, then
         * add the groups contribution from the multi-group fission source
         */
        virtual void fission( const ArrayF& fs, int ig );

        /**
         * Add the contribution from in-scattering from other groups. At some
         * point, ill play with upscattering iterations, but for now KISS.
         */
        virtual void in_scatter( size_t ig );

        /**
         * \brief Add a one-group auxiliary source
         * This adds some arbitrary source to the current group. Bear in mind
         * that the source definitiion starts with the MG fission source, then
         * contributions get tacked on from there.
         */
        void auxiliary( ArrayF aux ) {
            assert( source_1g_.size() == aux.size() );
            source_1g_ += aux;
        }

        /**
         * \brief Add self-scatter source
         * Add a contribution due to self-scatter within the current group,
         * returning the final source. This is usually called several times by a
         * sweeper in its "inner" iterations, and therefore does not mutate the
         * interal representation of the source, but instead returns the result
         * to the caller through the qbar argument.
         */
        virtual void self_scatter( size_t ig, ArrayF& flux_1g, 
                ArrayF& qbar ) const;

        /**
         * \brief Return the number of regions for which the Source is defined.
         */
        size_t n_reg() const {
            return n_reg_;
        }

    protected:
        const XSMesh *xs_mesh_;
        size_t n_group_;

        // This is true if an external source has been specified. For now it's
        // initialized false.
        bool has_external_;

        // Single-group source
        ArrayF source_1g_;

        // Reference to the MG flux variable. Need this to do scattering
        // contributions, etc.
        const ArrayF& flux_;

        size_t n_reg_;
    };

    typedef std::shared_ptr<Source> SP_Source_t;
    typedef std::unique_ptr<Source> UP_Source_t;
}
