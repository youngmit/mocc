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

#include <iosfwd>
#include <memory>
#include <vector>
#include "util/blitz_typedefs.hpp"
#include "util/error.hpp"
#include "util/global_config.hpp"
#include "util/range.hpp"
#include "core/angular_quadrature.hpp"
#include "core/coarse_data.hpp"
#include "core/eigen_interface.hpp"
#include "core/output_interface.hpp"
#include "core/source.hpp"
#include "core/source_factory.hpp"
#include "core/source_isotropic.hpp"
#include "core/xs_mesh.hpp"
#include "core/xs_mesh_homogenized.hpp"

namespace mocc {
/**
 * \todo clean up these constructors. Would be nice for the (input, mesh)
 * version to be able to call the (input) version.
 */
class TransportSweeper : public HasOutput {
public:
    TransportSweeper(const pugi::xml_node &input, const CoreMesh &mesh,
                     MeshTreatment treatment);

    TransportSweeper(const pugi::xml_node &input);

    virtual ~TransportSweeper()
    {
    }

    /**
     * \brief Perform a transport sweep of the passed group.
     */
    virtual void sweep(int group) = 0;

    /**
     * \brief Initialize the solution variables (scalar, boundary flux,
     * etc.) to reasonable initial guesses.
     */
    virtual void initialize() = 0;

    /**
     * \brief Update the incoming boundary flux values.
     *
     * This alters the incoming angular flux values to reflect the state of
     * the associated \ref CoarseData.
     */
    virtual void update_incoming_flux() = 0;

    /**
     * \brief Return a const reference to the Sweeper's \ref
     * AngularQuadrature
     */
    const AngularQuadrature &ang_quad() const
    {
        return ang_quad_;
    }

    /**
     * \brief Return a vector containing the pin-homogenized multi-group
     * scalar flux. The values in the vector are ordered group-major.
     */
    ArrayB2 get_pin_flux() const;

    /**
     * \brief Produce pin-homogenized scalar flux for the specified group
     * and store in the passed array.
     */
    virtual void get_pin_flux_1g(int ig, ArrayB1 &flux) const = 0;

    /**
     * \brief Project a single-group pin mesh-homogenized flux to the fine
     * mesh. Return the residual.
     */
    virtual real_t set_pin_flux_1g(int group, const ArrayB1 &pin_flux) = 0;

    /**
     * \brief Project a multi-group pin mesh-homogenized flux to the fine
     * mesh. Return the residual.
     */
    real_t set_pin_flux(const ArrayB2 &pin_flux)
    {
        real_t e = 0.0;
        for (int ig = 0; ig < (int)n_group_; ig++) {
            ArrayB1 flux_1g = pin_flux(blitz::Range::all(), ig);
            real_t e_g      = this->set_pin_flux_1g(ig, flux_1g);
            e += e_g * e_g;
        }
        return std::sqrt(e);
    }

    /**
     * \brief Return a const reference to the multi-group flux
     */
    const ArrayB2 &flux() const
    {
        return flux_;
    }

    /**
     * \brief Return a reference to the multi-group flux
     */
    ArrayB2 &flux()
    {
        return flux_;
    }

    /**
     * \brief Given the current estimate of a system eigenvalue, calculate
     * the group-independent fission source and store in the passed array
     */
    virtual void calc_fission_source(real_t k, ArrayB1 &fission_source) const;

    /**
     * \brief Construct and return a source object which conforms to the
     * sweeper.
     *
     * For now, default to the isotropic MoC Source type, \ref
     * SourceIsotropic.
     */
    virtual UP_Source_t create_source(const pugi::xml_node &input) const
    {
        try {
            auto source =
                SourceFactory(input, n_reg_, xs_mesh_.get(), this->flux());

            return source;
        } catch (Exception e) {
            std::cerr << e.what() << std::endl;
            throw EXCEPT("Failed to create source.");
        }
    }

    /**
     * \brief Return a shared pointer to a homogenized XS mesh.
     *
     * This is polymorphic, because some sweepers already operate on a
     * homogenized mesh, and there is no need to generate a new one.
     */
    virtual SP_XSMeshHomogenized_t get_homogenized_xsmesh() = 0;

    virtual int n_reg() const
    {
        return n_reg_;
    }

    /**
     * \brief Return the number of energy groups
     */
    int n_group() const
    {
        return n_group_;
    }

    /**
     * \brief Return a reference to the sweeper's XSMesh
     */
    const XSMesh &xs_mesh() const
    {
        return *(xs_mesh_.get());
    }

    /**
     * \brief Return a shared pointer to the sweeper's XSMesh. Use with
     * caution
     */
    SP_XSMesh_t get_xs_mesh()
    {
        return xs_mesh_;
    }

    /**
     * \brief Return a reference to the CoreMesh
     */
    const CoreMesh &mesh() const
    {
        return *core_mesh_;
    }

    /**
     * \brief Subscript and return a specific flux value
     */
    real_t flux(int ig, int ireg) const
    {
        assert(ig < (int)n_group_);
        assert(ireg < (int)n_reg_);

        return flux_(ireg, ig);
    }

    /**
     * \brief Assign a CoarseData object to the sweeper, allowing it to
     * store currents and such.
     */
    virtual void set_coarse_data(CoarseData *cd)
    {
        assert(cd);
        coarse_data_ = cd;
    }

    /**
     * \brief Associate the sweeper with a source.
     *
     * This is usually done by something like the \ref FixedSourceSolver.
     */
    virtual void assign_source(Source *source)
    {
        assert(source != nullptr);
        source_ = source;
    }

    /**
     * \brief Store the current flux as the old flux
     */
    virtual void store_old_flux()
    {
        flux_old_ = flux_;
        return;
    }

    /**
     * \brief Compute a flux residual between the current state of the flux
     * and the old flux. Defaults to L-2 norm.
     */
    virtual real_t flux_residual() const;

    /**
     * \brief Compute the total fission source based on the current or previous
     * state of the flux
     *
     * \param old whether to use the previous-iteration flux. Default: \c false
     */
    virtual real_t total_fission(bool old = false) const;

    /**
     * \brief Return a 3-D array containing normalized pin powers
     *
     * The nature of the normalization is somewhat up in the air. There are
     * different ways to do this. For instance, in the case of non-uniform
     * plane thicknesses and uniform power distribution, should a thicker
     * plane have a greater normalized pin power than a shorter plane? If
     * the pin volumes are not uniform (e.g. annular fuel), should a smaller
     * pin have less normalized power. Generally we would say "no" to the
     * former and "yes" to the latter, but this is pretty arbitrary, so...
     *
     * \note for now the normalization used is really simple! All values are
     * normalized uniformly such that the sum of all powers equals the
     * number of elements in the array.
     */
    virtual ArrayB3 pin_powers() const;

    /**
     * \brief Return a const reference to the region volumes
     */
    const VecF &volumes() const
    {
        return vol_;
    }

protected:
    const CoreMesh *core_mesh_;

    SP_XSMesh_t xs_mesh_;

    int n_reg_;
    int n_group_;
    std::vector<int> groups_;

    Source *source_;

    // Multi-group scalar flux
    ArrayB2 flux_;

    // Previous value of the MG scalar flux
    ArrayB2 flux_old_;

    // Region volumes
    VecF vol_;

    AngularQuadrature ang_quad_;

    // Reference to the CoarseData object that should be used to store
    // coarse mesh values. This is passed in from above.
    CoarseData *coarse_data_;

    // Total number of calls to sweep in the lifetime of the sweeper. Should
    // be n_group_ time the number of outer iterations
    int n_sweep_;

    // Total number of inner iteration sweeps
    int n_sweep_inner_;

    // Do incoming flux updates?
    bool do_incoming_update_;
};

typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
