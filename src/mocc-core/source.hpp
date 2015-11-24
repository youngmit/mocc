#pragma once

#include <iostream>

#include <blitz/array.h>

#include "blitz_typedefs.hpp"
#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    class Source {
    public:
        Source( int nreg, const XSMesh* xs_mesh, const ArrayB2& flux );

        virtual ~Source(){
        }

        /**
         * \brief Initialize the group source.
         *
         * If there is an external source specified, initialize to that.
         * Otherwise, initialize the external source.
         */
        virtual void initialize_group( int ig );

        /**
         * \brief Add the group's contribution from the multi-group fission
         * source.
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
        void auxiliary( const ArrayB1 &aux );

        /**
         * \brief Add self-scatter source
         * Add a contribution due to self-scatter within the current group,
         * returning the final source. This is usually called several times by a
         * sweeper in its "inner" iterations, and therefore does not mutate the
         * interal representation of the source, but instead returns the result
         * to the caller through the qbar argument.
         *
         * Have to pass flux_1g in for now, because it is assumed to be updated
         * in the sweeper inner iterations outside of the MG flux array. At some
         * point, it might be prudent to update it in-place, at which point,
         * this flux_1g argument here would be superfluous.
         */
        virtual void self_scatter( size_t ig, const ArrayB1& flux_1g,
                ArrayF& qbar ) const;

        /**
         * \brief Return the number of regions for which the Source is defined.
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * \brief Add an external source from an XML node.
         */
        void add_external( const pugi::xml_node &input );

        /**
         * \brief Return the source for the indexed cell for the current state of the
         * source.
         *
         * This will never include self-scatter. It will only have contributions
         * that have been added thus far.
         */
        real_t operator[]( size_t i ) const {
            return source_1g_[i];
        }

        /**
         * \copydoc Source::operator[]
         */
        real_t& operator[]( size_t i ) {
            return source_1g_[i];
        }


        friend std::ostream& operator<<(std::ostream &os, const Source &src) {
            std::cout << src.source_1g_ << std::endl;
            return os;
        }

    protected:
        const XSMesh *xs_mesh_;
        size_t n_group_;

        // This is true if an external source has been specified. For now it's
        // initialized false.
        bool has_external_;

        // The external source, if set
        ArrayB2 external_source_;

        // Single-group source. We use the Eigen storage class so that it can be
        // used directly as a source vector in a linear system.
        VectorX source_1g_;

        // Reference to the MG flux variable. Need this to do scattering
        // contributions, etc.
        const ArrayB2& flux_;

        size_t n_reg_;
    };

    typedef std::shared_ptr<Source> SP_Source_t;
    typedef std::unique_ptr<Source> UP_Source_t;
}
