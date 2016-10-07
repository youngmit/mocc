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
#include "sn_sweeper_factory_cdd.hpp"

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
    "dump_corrections",
    "update_incoming",
    "cycle"};
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
      prev_moc_flux_(sn_sweeper_->n_group(),
                     mesh_.n_reg(MeshTreatment::PIN_PLANE)),
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

    tl_ = 0.0;

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

    // Do an early Sn sweep if V-cycle is enabled
    if (v_cycle_) {
        sn_sweeper_->sweep(group);

        ArrayB1 sn_flux(mesh_.n_reg(MeshTreatment::PIN_PLANE));
        sn_sweeper_->get_pin_flux_1g(group, sn_flux, MeshTreatment::PIN_PLANE);

        // Check for negative fluxes on the Sn mesh
        int n_neg = 0;
        for (auto &v : sn_flux) {
            if (v < 0.0) {
                n_neg++;
                v = 0.0;
            }
        }
        if (n_neg > 0) {
            LogScreen << "Corrected " << n_neg
                      << " negative fluxes in Sn projection"
                      << "\n";
        }
        moc_sweeper_.set_pin_flux_1g(group, sn_flux, MeshTreatment::PIN_PLANE);
    }

    // MoC Sweeper
    bool do_moc =
        ((i_outer_ + 1) > n_inactive_moc_) && ((i_outer_ % moc_modulo_) == 0);
    if (do_moc) {
        moc_sweeper_.sweep(group);

        int n_negative  = 0;
        int n_NaN       = 0;
        const auto flux = moc_sweeper_.flux()(blitz::Range::all(), group);
        for (const auto &v : flux) {
            if (v < 0.0) {
                n_negative++;
            }
            if (v != v) {
                n_NaN++;
            }
        }
        if (n_negative > 0) {
            LogScreen << n_negative << " negative MoC fluxes in group " << group
                      << "\n";
        }
        if (n_NaN > 0) {
            LogScreen << n_NaN << " NaN MoC fluxes in group " << group << "\n";
        }
    }

    ArrayB1 prev_moc_flux = prev_moc_flux_(group, blitz::Range::all());
    moc_sweeper_.get_pin_flux_1g(group, prev_moc_flux,
                                 MeshTreatment::PIN_PLANE);

    if (do_mocproject_) {
        sn_sweeper_->set_pin_flux_1g(group, prev_moc_flux);
    }

    // Sn sweeper
    sn_sweeper_->sweep(group);

    ArrayB1 sn_flux(mesh_.n_reg(MeshTreatment::PIN_PLANE));
    sn_sweeper_->get_pin_flux_1g(group, sn_flux, MeshTreatment::PIN_PLANE);

    if (do_snproject_) {
        // Check for negative fluxes on the Sn mesh
        int n_neg = 0;
        for (auto &v : sn_flux) {
            if (v < 0.0) {
                n_neg++;
                v = 0.0;
            }
        }
        if (n_neg > 0) {
            LogScreen << "Corrected " << n_neg
                      << " negative fluxes in Sn projection"
                      << "\n";
        }
        moc_sweeper_.set_pin_flux_1g(group, sn_flux, MeshTreatment::PIN_PLANE);
    }

    // Compute Sn-MoC residual
    real_t residual = 0.0;
    for (int i = 0; i < (int)prev_moc_flux.size(); i++) {
        real_t diff = prev_moc_flux(i) - sn_flux(i);
        residual += diff * diff;
        sn_resid_(group, i) = diff;
    }
    residual = sqrt(residual) / mesh_.n_pin();

    LogScreen << "MoC/Sn residual: " << residual;
    if (sn_resid_norm_[group].size() > 0) {
        LogScreen << "   \t" << residual / sn_resid_norm_[group].back();
    }
    LogScreen << "\n";

    sn_resid_norm_[group].push_back(residual);
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::initialize()
{
    sn_sweeper_->initialize();
    moc_sweeper_.initialize();
}

////////////////////////////////////////////////////////////////////////////////
void PlaneSweeper_2D3D::get_pin_flux_1g(int ig, ArrayB1 &flux,
                                        MeshTreatment treatment) const
{
    // for now, only support PIN. We should only be calling this from places
    // above the sweeper itself, such as Eigen solver or in CMFD, so shouldnt
    // need the others.
    assert(treatment == MeshTreatment::PIN);

    if (expose_sn_) {
        sn_sweeper_->get_pin_flux_1g(ig, flux, MeshTreatment::PIN);
    } else {
        moc_sweeper_.get_pin_flux_1g(ig, flux, MeshTreatment::PIN);
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

    try {
        file.create_link("/Sn/xsmesh", "/xsmesh");
    } catch (...) {
        throw EXCEPT("Failed to create Sn xsmesh link");
    }
    try {
        file.create_link("/Sn/ang_quad", "/ang_quad");
    } catch (...) {
        throw EXCEPT("Failed to create Sn ang_quad link");
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
            h5g.write(setname.str(), flux(ig, blitz::Range::all()), dims_moc);
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
    v_cycle_                = false;

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

    if (!input.attribute("cycle").empty()) {
        std::string cycle = input.attribute("cycle").value();
        if (cycle == "v") {
            v_cycle_ = true;
        } else if ((cycle == "saw") || (cycle == "sawtooth")) {
            v_cycle_ = false;
        } else {
            throw EXCEPT("Unrecognized cycle attribute");
        }
    }

    // Make sure that sn project is on if we are exposing sn
    if (expose_sn_ && !do_snproject_) {
        Warn(
            "Exposing Sn as global solver and not projecting to MoC. This will "
            "cause weirdness in the fission source normalization.");
    }

    // Throw a warning if TL is disabled
    if (!do_tl_) {
        Warn(
            "Transverse leakage is disabled. Are you sure that's what you "
            "want?");
    }

    LogFile << "2D3D Sweeper options:"
            << "\n";
    LogFile << "    Sn Projection: " << do_snproject_ << "\n";
    LogFile << "    MoC Projection: " << do_mocproject_ << "\n";
    LogFile << "    Expose Sn pin flux: " << expose_sn_ << "\n";
    LogFile << "    Keep original Sn quadrature: " << keep_sn_quad_ << "\n";
    LogFile << "    Transverse Leakage: " << do_tl_ << "\n";
    LogFile << "    Relaxation factor: " << relax_ << "\n";
    LogFile << "    Inactive MoC Outer Iterations: " << n_inactive_moc_ << "\n";
    LogFile << "    MoC sweep modulo: " << moc_modulo_ << "\n";
    LogFile << "    Apply Sn-MoC flux residual to CMFD updates: "
            << discrepant_flux_update_ << "\n";
    LogFile << "    Sweep cycle: ";
    if (v_cycle_) {
        LogFile << "V"
                << "\n";
    } else {
        LogFile << "Sawtooth"
                << "\n";
    }
}
}
} // Namespace mocc::cmdo
