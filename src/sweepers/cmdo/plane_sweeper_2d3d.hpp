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
#include "util/global_config.hpp"
#include "util/pugifwd.hpp"
#include "core/angular_quadrature.hpp"
#include "core/output_interface.hpp"
#include "sn/sn_sweeper_variant.hpp"
#include "correction_data.hpp"
#include "moc_sweeper_2d3d.hpp"
#include "sn_sweeper_cdd.hpp"
#include "sn_sweeper_factory_cdd.hpp"
#include "source_2d3d.hpp"

namespace mocc {
namespace cmdo {
/**
 * This is an implementation of the 2D3D method. Each plane is treated with
 * a 2-D MoC sweeper, which produces the correction factors needed to treat
 * the entire system with a 3-D corrected diamond difference Sn sweeper.
 */
class PlaneSweeper_2D3D : public TransportSweeper {
public:
    PlaneSweeper_2D3D(const pugi::xml_node &input, const CoreMesh &mesh);

    void sweep(int group) override final;

    void initialize() override final;

    /**
     * \brief \copybrief mocc::TransportSweeper::update_incoming_flux()
     *
     * This delegates to both contained sweepers.
     */
    void update_incoming_flux() override final
    {
        Warn("Incoming flux updates are not properly supported yet for 2D3D");

        // moc_sweeper_.update_incoming_flux();
        // sn_sweeper_->update_incoming_flux();

        return;
    }

    /**
     * \brief \copybrief TransportSweeper::get_pin_flux_1g()
     */
    void get_pin_flux_1g(int ig, ArrayB1 &flux,
                         MeshTreatment treatment) const override final;

    /**
     * \brief \copybrief TransportSweeper::set_pin_flux_1g()
     *
     * Delegate to the subbordinate \ref sn::SnSweeper and \ref
     * moc::MoCSweeper.  Return the error from the MoC sweeper.
     */
    real_t
    set_pin_flux_1g(int group, const ArrayB1 &pin_flux,
                    MeshTreatment treatment = MeshTreatment::PIN) override final
    {
        assert((treatment == MeshTreatment::PIN) ||
               (treatment == MeshTreatment::PIN_PLANE));
        sn_sweeper_->set_pin_flux_1g(group, pin_flux, treatment);

        real_t diff = 0.0;
        if (discrepant_flux_update_) {
            assert(false);
            ArrayB1 moc_pin_flux(pin_flux.size());
            for (int i = 0; i < (int)moc_pin_flux.size(); i++) {
                moc_pin_flux(i) =
                    std::max(0.0, pin_flux(i) - sn_resid_(group, i));
            }
        } else {
            diff = moc_sweeper_.set_pin_flux_1g(group, pin_flux,
                                                treatment);
        }

        return diff;
    }

    /**
     * \brief \copybrief HasOutput::output()
     */
    void output(H5Node &file) const override final;

    /**
     * \brief \copybrief TransportSweeper::assign_source()
     *
     * Associate the sweeper with a source. This has to do a little extra
     * work, since the Sn sweeper needs its own source.
     */
    virtual void assign_source(Source *source) override final
    {
        assert(source != nullptr);
        source_ = source;
        moc_sweeper_.assign_source(source);
        /// \todo this static_cast is scary. Maybe think about relaxing the
        /// ownership of the source by FSS and allow the TS to figure out
        /// the types more explicitly...
        Source_2D3D *s = static_cast<Source_2D3D *>(source);
        sn_sweeper_->assign_source(s->get_sn_source());
    }

    /**
     * \brief \copybrief TransportSweeper::create_source()
     *
     * Create a Source_2D3D object instead of the standard Source class.
     *
     * \todo This doesnt play nice with the existing source factory. This is
     * mostly because the source factory does a little too much setup
     * based on the input to the \<source /\> tag, which is better done by
     * the Source constructor. Clean this up
     */
    UP_Source_t create_source(const pugi::xml_node &input) const override final
    {
        std::cout << "creating 2d3d source" << std::endl;

        auto source = UP_Source_t(new Source_2D3D(moc_sweeper_, *sn_sweeper_));
        return source;
    }

    SP_XSMeshHomogenized_t get_homogenized_xsmesh() override final
    {
        return sn_sweeper_->get_homogenized_xsmesh();
    }

    int n_reg() const override final
    {
        if (expose_sn_) {
            return sn_sweeper_->n_reg();
        } else {
            return moc_sweeper_.n_reg();
        }
    }

    /**
     * \brief \copybrief TransportSweeper::n_reg_fission()
     *
     * This returns the sum of the sizes of the Sn and MoC sweeper. Calls to
     * \ref PlaneSweeper_2D3D::calc_fission_source() will store the respective
     * fission sources for the two sweepers one after the other.
     */
    int n_reg_fission() const override final
    {
        return sn_sweeper_->n_reg() + moc_sweeper_.n_reg();
    }

    /**
     * \brief \copybrief TransportSweeper::calc_fission_source()
     *
     * Override the default \ref TransportSweeper implementation to call the
     * method on one of the sub-sweepers.
     */
    void calc_fission_source(real_t k,
                             ArrayB1 &fission_source) const override final
    {
        assert((int)fission_source.size() ==
               moc_sweeper_.n_reg() + sn_sweeper_->n_reg());

        // Alias the appropriate regions of the passed fission source
        ArrayB1 sn_fission_source(
            fission_source(blitz::Range(0, sn_sweeper_->n_reg() - 1)));
        ArrayB1 moc_fission_source(
            fission_source(blitz::Range(sn_sweeper_->n_reg(), blitz::toEnd)));
        sn_sweeper_->calc_fission_source(k, sn_fission_source);
        moc_sweeper_.calc_fission_source(k, moc_fission_source);

        return;
    }

    /**
     * \brief \copybrief TransportSweeper::total_fission()
     *
     * Override the default \ref TransportSweeper implementation to call the
     * method on one of the sub-sweepers. For now, using the MoC
     * implementation, since it's the finer mesh, generally speaking
     */
    real_t total_fission(bool old) const override final
    {
        if (expose_sn_) {
            return sn_sweeper_->total_fission(old);
        } else {
            return moc_sweeper_.total_fission(old);
        }
    }

    /**
     * \brief \copybrief TransportSweeper::store_old_flux()
     *
     * Defer to the MoC and Sn sweepers.
     */
    void store_old_flux() override final
    {
        moc_sweeper_.store_old_flux();
        sn_sweeper_->store_old_flux();
        return;
    }

    /**
     * \brief \copybrief TransportSweeper::set_coarse_data()
     *
     * Delegate to subordinate sweepers.
     */
    void set_coarse_data(CoarseData *cd) override final
    {
        coarse_data_ = cd;
        moc_sweeper_.set_coarse_data(cd);
        sn_sweeper_->set_coarse_data(cd);
    }

private:
    // Secret constructor, delegated to by the public constructor. This only
    // exists so that the public ctor can store the result of the CDD factory in
    // a temporary to be used in this ctor's initializer list.
    PlaneSweeper_2D3D(const pugi::xml_node &input, const CoreMesh &mesh,
                      CDDPair_t cdd_pair);
    // Parse the various options from the XML
    void parse_options(const pugi::xml_node &input);
    // Calculate transverse leakage based on the state of the coarse_data_
    // and apply to the MoC sweeper's source.
    void add_tl(int group);

    const CoreMesh &mesh_;

    // The number of pins to consider in the MoC sweeper
    int n_pin_moc_;

    UP_SnSweeper_t sn_sweeper_;
    std::shared_ptr<CorrectionData> corrections_;
    MoCSweeper_2D3D moc_sweeper_;
    AngularQuadrature ang_quad_;
    // Pin-level transverse leakage. This is in coarse mesh ordering (at the
    // macroplane axial refinement) for easy output to HDF5, since it is
    // interacted with from the MoC side in expanded FSR space anyways
    ArrayB2 tl_;

    // L-2 norm of the Sn-MoC residuals by group sweep
    std::vector<VecF> sn_resid_norm_;

    // Actual Sn-MoC residuals
    ArrayB2 sn_resid_;

    // Pre-Sn projection moc flux. Useful for keeping track of MoC-Sn
    // residual
    ArrayB2 prev_moc_flux_;

    // Outer iteration index. Starts at -1 and is incremented whenever group
    // 0 is swept. This is kind of brittle.
    int i_outer_;

    // Options! Buttons and knobs!!!
    // When asking the sweeper for pin flux, which one?
    bool expose_sn_;
    // Project Sn flux to fine mesh after Sn sweeps? This is important when
    // no CMFD acceleration is used.
    bool do_snproject_;
    // Project MoC flux to Sn after MoC inners?
    bool do_mocproject_;
    // Whether we should replace the sn angular quadrature with the
    // modularized quadrature from MoC. Technically we need to do this
    // always, but turning this off sometimes is useful for debugging
    bool keep_sn_quad_;
    // Enable transverse leakage? In most cases this will prevent
    // convergence.
    bool do_tl_;
    // Number of outer iterations to skip MoC. Super experimental
    int n_inactive_moc_;
    int moc_modulo_;
    // Relaxation factor for the flux updates
    real_t relax_;
    // Whether to incorporate MoC/Sn error in CMFD flux update
    bool discrepant_flux_update_;
    // Write correction factors to HDF5 file?
    bool dump_corrections_;
    // Whether to use a sawtooth or V cycle in the sweep
    bool v_cycle_;
};
}
} // Namespace mocc::cmdo
