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

#include <blitz/array.h>
#include <iosfwd>
#include "util/blitz_typedefs.hpp"
#include "util/global_config.hpp"
#include "core/eigen_interface.hpp"
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
 * class). How the implementation will work is somewhat up in the air. There
 * is also a question as to how arbitrary sources will work with the
 * fixed-source solver, so for now the split between the base \ref Source
 * class and its derived \ref SourceIsotropic is a little arbitrary, and
 * should become more clear as Pn and arbitrary sources get implemented.
 */
class Source {
public:
    Source(int nreg, const XSMesh *xs_mesh, const ArrayB2 &flux);

    virtual ~Source()
    {
    }

    /**
     * \brief Initialize the group source.
     *
     * If there is an external source specified, initialize to that.
     * Otherwise, initialize the external source.
     */
    virtual void initialize_group(int ig);

    /**
     * \brief Add the group's contribution from the multi-group fission
     * source.
     */
    virtual void fission(const ArrayB1 &fs, int ig);

    /**
     * Add the contribution from in-scattering from other groups. At some
     * point, ill play with up-scattering iterations, but for now KISS.
     */
    virtual void in_scatter(size_t ig);

    /**
     * \brief Add a one-group auxiliary source
     *
     * This adds some arbitrary source to the current group. Bear in mind
     * that the source definition starts with the MG fission source, then
     * contributions get tacked on from there.
     */
    void auxiliary(const ArrayB1 &aux);

    /**
     * \brief Add self-scatter source
     *
     * Add a contribution due to self-scatter within the current group,
     * producing the final source. This is usually called several times by a
     * sweeper in its "inner" iterations, and therefore does not mutate the
     * internal representation of the source stored in source_1g_, but
     * instead stores the result in a different location, which should be
     * accessible via the \ref get_transport() method.
     */
    virtual void self_scatter(size_t ig, const ArrayB1 &xstr = ArrayB1(0)) = 0;

    /**
     * \brief Add self-scatter source
     *
     * This is different from the above self-scatter function in that we only
     * want to output the 1g source with the self-scattering source for use in
     * MMS.
     */

    virtual void self_scatter_for_MMS(size_t ig, const ArrayB1 &xstr = ArrayB1(0)) = 0;


    /**
     * \brief Return the number of regions for which the Source is defined.
     */
    int n_reg() const
    {
        return n_reg_;
    }

    /**
     * \brief Add an external source from an XML node.
     */
    void add_external(const pugi::xml_node &input);

    /**
     * \brief Scale the source by some weighting values.
     *
     * Use this if the client code desires a total source [n] for each
     * region, rather than a specific source [nn/cm-3]. Make sure to use
     * this right before using the source, after adding all contributions.
     */
    void scale(const VecF &v)
    {
        assert((int)v.size() == n_reg_);
        for (int i = 0; i < n_reg_; i++) {
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
    real_t operator[](size_t i) const
    {
        return source_1g_[i];
    }

    /**
     * \copydoc Source::operator[]
     */
    real_t &operator[](size_t i)
    {
        return source_1g_[i];
    }

    /**
     * \brief Return a const reference to the actual 1G source.
     */
    const VectorX &get() const
    {
        return source_1g_;
    }

    /**
     * \brief Return a const reference to the 1G source with self-scattering source.
     */

    const VectorX &get_source_1g_with_self_scat(int iang) const
    {
        return source_1g_with_self_scat_;
    }

    /**
     * \brief Return a const reference to the source, as it should be used
     * in a transport sweeper.
     *
     * \param [in] iang the angle index to get a source for. Not used for
     * isotropic sources, but necessary for angle-dependent sources.
     */
    virtual const VectorX &get_transport(int iang) const = 0;

    // Get group wise fission source, the purpose for now is for output the
    // Multigroup fission source
    const ArrayB2 get_mg_fission_source(const ArrayB1 &fs);

    friend std::ostream &operator<<(std::ostream &os, const Source &src);

protected:
    /**
     * This struct stores the current state of the \ref Source. While not
     * absolutely necessary, it keeps client code from doing stupid things.
     */
    struct State {
        bool has_fission : 1;
        bool has_inscatter : 1;
        bool is_scaled : 1;
        void reset()
        {
            has_fission   = false;
            has_inscatter = false;
            is_scaled     = false;
        }
    };

    State state_;

    /**
     * Reference to a compatible \ref TransportSweeper \ref XSMesh.
     */
    const XSMesh *xs_mesh_;

    int n_group_;
    int n_reg_;

    // This is true if an external source has been specified. For now it's
    // initialized false.
    bool has_external_;

    // The external source, if set
    ArrayB2 external_source_;

    // Single-group source. We use the Eigen storage class so that it can be
    // used directly as a source vector in a linear system.
    VectorX source_1g_;

    // The above source_1g_ does not include self-scattering source.
    // The following one does, mainly for use of MMS
    VectorX source_1g_with_self_scat_;

    // Reference to the MG flux variable. Need this to do scattering
    // contributions, etc.
    const ArrayB2 &flux_;
};

typedef std::shared_ptr<Source> SP_Source_t;
typedef std::unique_ptr<Source> UP_Source_t;
}
