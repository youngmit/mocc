/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include <iostream>

#include <blitz/array.h>

#include "blitz_typedefs.hpp"
#include "eigen_interface.hpp"
#include "global_config.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    /**
     * The \ref Source class provides functionality for defining a particle
     * source, accumulated from commonly-used origins (in-scattering, fission,
     * self-scatter, etc.). A \ref Source object can be thought of as a simple
     * state machine, in which source contributions for a single group (and
     * potentially, angle) are added one-at-a-time. The "state" is reset with a
     * call to \ref Source::initialize_group(), which clears all contributions
     * to the internally-stored single-group source by initializing it to a
     * prescribed external source (if one is defined), or to zero. Following
     * group initialization, the other methods can be invoked by the client code
     * to add other contributions as needed. When all contributions have been
     * added, the client can then request a source for the current group and
     * desired angle.
     *
     * \note At some point, it may be nice to incorporate Pn scattering. Most of
     * the functionality for Pn scattering would end up in the Source class
     * hierarchy (probably by introducing a Pn variant of the \ref Source base
     * class). How the implementation will occur is somewhat up in the air.
     * Simplest approach would be to add a new virtual method (which just
     * <tt>return</tt>s on non-Pn sources), which accepts angular flux as input
     * somehow.  This could be done by passing in an FSR-dependent angular flux
     * after each angle sweep and have the \ref Source contribute to angular
     * flux moments internally. This probably isn't the most computationally
     * efficient method, though. Ultimately, flux moments need to be stored
     * somewhere, and the \ref Source is a great candidate; perhaps exposing a
     * reference to those moments to the \ref TransportSweeper and coming up
     * with a slick way for the \ref TransportSweeper to interact with those
     * moments in an as-needed way would work better. Cross that river when we
     * get to it, I suppose...
     */
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
        virtual void fission( const ArrayB1& fs, int ig );

        /**
         * Add the contribution from in-scattering from other groups. At some
         * point, ill play with up-scattering iterations, but for now KISS.
         */
        virtual void in_scatter( size_t ig );

        /**
         * \brief Add a one-group auxiliary source
         * This adds some arbitrary source to the current group. Bear in mind
         * that the source definition starts with the MG fission source, then
         * contributions get tacked on from there.
         */
        void auxiliary( const ArrayB1 &aux );

        /**
         * \brief Add self-scatter source
         * Add a contribution due to self-scatter within the current group,
         * returning the final source. This is usually called several times by a
         * sweeper in its "inner" iterations, and therefore does not mutate the
         * internal representation of the source, but instead returns the result
         * to the caller through the qbar argument.
         */
        virtual void self_scatter( size_t ig ) = 0;

        /**
         * \brief Return the number of regions for which the Source is defined.
         */
        size_t n_reg() const {
            return n_reg_;
        }

        /**
         * \brief Define whether this \ref Source should scale itself by the
         * transport cross section.
         */
        void set_scale_transport( bool scale ) {
            scale_transport_ = scale;
        }

        /**
         * \brief Add an external source from an XML node.
         */
        void add_external( const pugi::xml_node &input );

        /**
         * \brief Scale the source by some weighting values.
         *
         * Use this if the client code desires a total source [n] for each
         * region, rather than a specific source [nn/cm-3]. Make sure to use
         * this right before using the source, after adding all contributions.
         */
        void scale( const VecF &v ) {
            assert( v.size() == n_reg_ );
            for( int i=0; i<(int)n_reg_; i++ ) {
                source_1g_(i) *= v[i];
            }
        }

        /**
         * \brief Return the source for the indexed cell for the current state
         * of the source.
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

        /**
         * \brief Return a const reference to the actual 1G source.
         */
        const VectorX& get() const {
            return source_1g_;
        }

        /**
         * \brief Return a const reference to the source, as it should be used
         * in a transport sweeper.
         *
         * \param [in] iang the angle index to get a source for. Not used for
         * isotropic sources, but necessary for angle-dependent sources.
         */
        virtual const VectorX& get_transport( int iang ) const = 0;

        friend std::ostream& operator<<(std::ostream &os, const Source &src) {
            std::cout << src.source_1g_ << std::endl;
            return os;
        }

    protected:
        /**
         * This struct stores the current state of the \ref Source. While not
         * absolutely necessary, it keeps client code from doing stupid things.
         */
        struct State {
            bool has_fission : 1;
            bool has_inscatter : 1;
            bool is_scaled : 1;
            void reset() {
                has_fission = false;
                has_inscatter = false;
                is_scaled = false;
            }
        };

        State state_;

        /**
         * Reference to a compatible \ref TransportSweeper \ref XSMesh.
         */
        const XSMesh *xs_mesh_;

        size_t n_group_;

        // This is true if an external source has been specified. For now it's
        // initialized false.
        bool has_external_;

        // Should we scale the source by the transport cross section? This is a
        // useful optimization for MoC, but not necessarily other sweeper types.
        bool scale_transport_;

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
