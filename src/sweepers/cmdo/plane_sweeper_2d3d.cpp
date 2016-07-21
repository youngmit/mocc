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

#include "plane_sweeper_2d3d.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include "util/error.hpp"
#include "util/range.hpp"
#include "util/validate_input.hpp"

using mocc::sn::SnSweeper;

namespace {
const std::vector<std::string> recognized_attributes = {
    "type",
    "expose_sn",
    "sn_project",
    "moc_project",
    "tl",
    "inactive_moc",
    "moc_modulo",
    "preserve_sn_quadrature",
    "relax",
    "discrepant_flux_update",
    "dump_corrections"};
}

namespace mocc {
namespace cmdo {
////////////////////////////////////////////////////////////////////////////////
PlaneSweeper_2D3D::PlaneSweeper_2D3D(const pugi::xml_node &input,
                                     const CoreMesh &mesh)
    : PlaneSweeper_2D3D(input, mesh,
                        SnSweeperFactory_CDD(input.child("sn_sweeper"), mesh))
{
    return;
}
////////////////////////////////////////////////////////////////////////////////
PlaneSweeper_2D3D::PlaneSweeper_2D3D(const pugi::xml_node &input,
                                     const CoreMesh &mesh, CDDPair_t cdd_pair)
    : TransportSweeper(input),
      mesh_(mesh),
      n_pin_moc_(std::accumulate(
          mesh_.macroplanes().begin(), mesh_.macroplanes().end(), 0,
          [](int n, const MacroPlane p) { return n + p.plane->n_pin(); })),
      sn_sweeper_(std::move(cdd_pair.first)),
      corrections_(cdd_pair.second),
      moc_sweeper_(input.child("moc_sweeper"), mesh),
      ang_quad_(moc_sweeper_.get_ang_quad()),
      tl_(sn_sweeper_->n_group(), n_pin_moc_),
      sn_resid_norm_(sn_sweeper_->n_group()),
      sn_resid_(sn_sweeper_->n_group(), mesh_.n_pin()),
      prev_moc_flux_(sn_sweeper_->n_group(), mesh_.n_pin()),
      i_outer_(-1)
{
    validate_input(input, recognized_attributes);
    this->parse_options(input);
    core_mesh_ = &mesh;

    xs_mesh_ = moc_sweeper_.get_xs_mesh();
    flux_.reference(moc_sweeper_.flux());
    vol_ = moc_sweeper_.volumes();

    n_reg_   = moc_sweeper_.n_reg();
    n_group_ = xs_mesh_->n_group();
    groups_  = Range(n_group_);

    auto sn_xs_mesh = sn_sweeper_->get_homogenized_xsmesh();
    assert(corrections_);
    moc_sweeper_.set_coupling(corrections_, sn_xs_mesh,
                              sn_sweeper_->expanded_xs());

    if (!keep_sn_quad_) {
        sn_sweeper_->set_ang_quad(ang_quad_);
    }

    sn_sweeper_->get_homogenized_xsmesh()->set_flux(moc_sweeper_.flux());

    coarse_data_ = nullptr;

    return;
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::sweep(int group)
{
    if (!coarse_data_) {
        throw EXCEPT("CMFD must be enabled to do 2D3D.");
    }

    /// \todo do something less brittle
    if (group == 0) {
        i_outer_++;
    }

    // Calculate transverse leakage source
    if (do_tl_) {
        this->add_tl(group);
    }

    // MoC Sweeper
    bool do_moc =
        ((i_outer_ + 1) > n_inactive_moc_) && ((i_outer_ % moc_modulo_) == 0);
    if (do_moc) {
        moc_sweeper_.sweep(group);

        int n_negative  = 0;
        const auto flux = moc_sweeper_.flux()(blitz::Range::all(), group);
        for (const auto &v : flux) {
            if (v < 0.0) {
                n_negative++;
            }
        }
        if (n_negative > 0) {
            LogScreen << n_negative << " negative fluxes in group " << group
                      << std::endl;
        }
    }

    ArrayB1 prev_moc_flux = prev_moc_flux_(group, blitz::Range::all());
    moc_sweeper_.get_pin_flux_1g(group, prev_moc_flux);

    if (do_mocproject_) {
        sn_sweeper_->set_pin_flux_1g(group, prev_moc_flux);
    }

    // Sn sweeper
    sn_sweeper_->sweep(group);

    if (do_snproject_) {
        ArrayB1 sn_flux(mesh_.n_pin());
        sn_sweeper_->get_pin_flux_1g(group, sn_flux);
        moc_sweeper_.set_pin_flux_1g(group, sn_flux);
    }

    // Compute Sn-MoC residual
    real_t residual = 0.0;
    for (size_t i = 0; i < prev_moc_flux.size(); i++) {
        residual += (prev_moc_flux(i) - sn_sweeper_->flux(group, i)) *
                    (prev_moc_flux(i) - sn_sweeper_->flux(group, i));
        sn_resid_(group, (int)i) =
            sn_sweeper_->flux(group, i) - prev_moc_flux(i);
    }
    residual = sqrt(residual) / mesh_.n_pin();

    LogScreen << "MoC/Sn residual: " << residual;
    if (sn_resid_norm_[group].size() > 0) {
        LogScreen << "   \t" << residual / sn_resid_norm_[group].back();
    }
    LogScreen << std::endl;

    sn_resid_norm_[group].push_back(residual);
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::initialize()
{
    sn_sweeper_->initialize();
    moc_sweeper_.initialize();
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::get_pin_flux_1g(int ig, ArrayB1 &flux) const
{
    if (expose_sn_) {
        sn_sweeper_->get_pin_flux_1g(ig, flux);
    } else {
        moc_sweeper_.get_pin_flux_1g(ig, flux);
    }
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::add_tl(int group)
{
    assert(coarse_data_);
    ArrayB1 tl_fsr(n_reg_);

    blitz::Array<real_t, 1> tl_g = tl_(group, blitz::Range::all());

    int iplane   = 0;
    int ireg_pin = 0;
    int ipin     = 0;
    for (const auto &mplane : mesh_.macroplanes()) {
        real_t dz = mplane.height;
        for (const auto &mpin : mplane) {
            Position pos  = mesh_.pin_position(ipin);
            pos.z         = mplane.iz_min;
            size_t icell  = mesh_.coarse_cell(pos);
            int surf_down = mesh_.coarse_surf(icell, Surface::BOTTOM);
            pos.z         = mplane.iz_max;
            icell         = mesh_.coarse_cell(pos);
            int surf_up   = mesh_.coarse_surf(icell, Surface::TOP);
            // Get index for storing in tl_ array
            pos.z       = iplane;
            int icoarse = mesh_.coarse_cell(pos);

            real_t j_up   = coarse_data_->current(surf_up, group);
            real_t j_down = coarse_data_->current(surf_down, group);
            tl_g(icoarse) =
                tl_g(icoarse) * (1.0 - relax_) + relax_ * (j_down - j_up) / dz;
            for (int ir = 0; ir < mpin->n_reg(); ir++) {
                tl_fsr(ireg_pin) = tl_g(icoarse);
                ireg_pin++;
            }
            ipin++;
        }
        iplane++;
    }

    // Hand the transverse leakage to the MoC sweeper.
    moc_sweeper_.apply_transverse_leakage(group, tl_fsr);
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::output(H5Node &file) const
{
    // Put the Sn data in its own location
    {
        auto g = file.create_group("/Sn");
        sn_sweeper_->output(g);
    }

    // Put the MoC data in its own location
    {
        auto g = file.create_group("/MoC");
        moc_sweeper_.output(g);
    }

    VecI dims_moc;
    dims_moc.push_back(mesh_.macroplanes().size());
    dims_moc.push_back(mesh_.ny());
    dims_moc.push_back(mesh_.nx());
    VecI dims;
    dims.push_back(mesh_.nz());
    dims.push_back(mesh_.ny());
    dims.push_back(mesh_.nx());

    // Write out the Sn-MoC residual convergence
    file.create_group("/SnResid");
    for (int g = 0; g < n_group_; g++) {
        std::stringstream setname;
        setname << "/SnResid/" << std::setfill('0') << std::setw(3) << g;
        VecI niter(1, sn_resid_norm_[g].size());
        file.write(setname.str(), sn_resid_norm_[g], niter);
    }

    {
        auto flux = prev_moc_flux_.copy();
        Normalize(flux.begin(), flux.end());
        auto h5g = file.create_group("moc_flux");
        for (const auto &ig : groups_) {
            std::stringstream setname;
            setname << std::setfill('0') << std::setw(3) << ig + 1;
            h5g.write(setname.str(), flux(ig, blitz::Range::all()), dims);
        }
    }

    // Write out the transverse leakages
    {
        auto group = file.create_group("/transverse_leakage");
        for (int g = 0; g < n_group_; g++) {
            std::stringstream setname;
            setname << std::setfill('0') << std::setw(3) << g;

            const auto tl_slice = tl_((int)g, blitz::Range::all());

            group.write(setname.str(), tl_slice.begin(), tl_slice.end(),
                        dims_moc);
        }
    }

    // Write out the correction factors
    if (dump_corrections_) {
        corrections_->output(file);
    }
}

////////////////////////////////////////////////////////////////////////////////
// At some point it might be nice to make the options const and initialized
// them in the initializer list, then just check for validity later. This if
// fine for now.
void PlaneSweeper_2D3D::parse_options(const pugi::xml_node &input)
{
    // Set defaults for everything
    expose_sn_              = false;
    do_snproject_           = false;
    do_mocproject_          = false;
    keep_sn_quad_           = false;
    do_tl_                  = true;
    n_inactive_moc_         = 0;
    moc_modulo_             = 1;
    relax_                  = 1.0;
    discrepant_flux_update_ = false;
    dump_corrections_       = false;

    // Override with entries in the input node
    if (!input.attribute("expose_sn").empty()) {
        expose_sn_ = input.attribute("expose_sn").as_bool();
    }
    if (!input.attribute("sn_project").empty()) {
        do_snproject_ = input.attribute("sn_project").as_bool();
    }
    if (!input.attribute("moc_project").empty()) {
        do_mocproject_ = input.attribute("moc_project").as_bool();
    }
    if (!input.attribute("tl").empty()) {
        do_tl_ = input.attribute("tl").as_bool();
    }
    if (!input.attribute("inactive_moc").empty()) {
        n_inactive_moc_ = input.attribute("inactive_moc").as_int();
    }
    if (!input.attribute("moc_modulo").empty()) {
        moc_modulo_ = input.attribute("moc_modulo").as_int();
    }
    if (!input.attribute("preserve_sn_quadrature").empty()) {
        keep_sn_quad_ = input.attribute("preserve_sn_quadrature").as_bool();
    }
    if (!input.attribute("relax").empty()) {
        relax_ = input.attribute("relax").as_double(1.0);
    }
    if (!input.attribute("discrepant_flux_update").empty()) {
        discrepant_flux_update_ =
            input.attribute("discrepant_flux_update").as_bool();
    }
    dump_corrections_ = input.attribute("dump_corrections").as_bool(false);

    // Throw a warning if TL is disabled
    if (!do_tl_) {
        Warn("Transverse leakage is disabled. Are you sure that's what you "
             "want?");
    }

    // Make sure that if we are doing expose_sn, we arent also trying to do
    // MoC. Doesnt work right now.
    if (expose_sn_) {
        // Cheat and peek into the MoC tag
        int n_inner = input.child("moc_sweeper").attribute("n_inner").as_int(0);
        if (n_inner > 0) {
            Warn("Probably shouldn't expose the Sn sweeper while "
                 "doing MoC sweeps");
        }
    }

    LogFile << "2D3D Sweeper options:" << std::endl;
    LogFile << "    Sn Projection: " << do_snproject_ << std::endl;
    LogFile << "    MoC Projection: " << do_mocproject_ << std::endl;
    LogFile << "    Expose Sn pin flux: " << expose_sn_ << std::endl;
    LogFile << "    Keep original Sn quadrature: " << keep_sn_quad_
            << std::endl;
    LogFile << "    Transverse Leakage: " << do_tl_ << std::endl;
    LogFile << "    Relaxation factor: " << relax_ << std::endl;
    LogFile << "    Inactive MoC Outer Iterations: " << n_inactive_moc_
            << std::endl;
    LogFile << "    MoC sweep modulo: " << moc_modulo_ << std::endl;
    LogFile << "    Apply Sn-MoC flux residual to CMFD updates: "
            << discrepant_flux_update_ << std::endl;
}
}
} // Namespace mocc::cmdo
