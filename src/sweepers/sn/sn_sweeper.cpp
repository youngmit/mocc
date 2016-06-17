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
#include "sn_sweeper.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include "pugixml.hpp"
#include "util/files.hpp"
#include "util/string_utils.hpp"

namespace {
using namespace mocc;
BC_Size_t boundary_helper(const Mesh &mesh)
{
    BC_Size_t bc_size = {(int)mesh.ny() * (int)mesh.nz(),
                         (int)mesh.nx() * (int)mesh.nz(),
                         (int)mesh.nx() * (int)mesh.ny()};
    return bc_size;
}
}

namespace mocc {
namespace sn {
SnSweeper::SnSweeper(const pugi::xml_node &input, const CoreMesh &mesh)
    : TransportSweeper(input),
      timer_(RootTimer.new_timer("Sn Sweeper", true)),
      timer_init_(timer_.new_timer("Initialization", true)),
      timer_sweep_(timer_.new_timer("Sweep", false)),
      timer_xsupdate_(timer_.new_timer("XS Update", false)),
      mesh_(mesh),
      bc_type_(mesh.boundary()),
      flux_1g_(),
      xstr_(mesh.n_pin()),
      bc_in_(mesh.mat_lib().n_group(), ang_quad_, bc_type_,
             boundary_helper(mesh)),
      bc_out_(1, ang_quad_, bc_type_, boundary_helper(mesh)),
      gs_boundary_(true)
{
    LogFile << "Constructing a base Sn sweeper" << std::endl;

    // Set up the cross-section mesh. If there is <data> specified, try to
    // use that, otherwise generate volume-weighted cross sections
    if (input.child("data").empty()) {
        xs_mesh_ = SP_XSMesh_t(new XSMeshHomogenized(mesh));
    }
    else {
        try {
            xs_mesh_ = SP_XSMesh_t(new XSMeshHomogenized(mesh, input));
        }
        catch (Exception e) {
            std::cerr << e.what() << std::endl;
            throw EXCEPT("Failed to create XSMesh for Sn Sweeper.");
        }
    }

    core_mesh_ = &mesh;
    n_reg_     = mesh.n_pin();
    n_group_   = xs_mesh_->n_group();
    flux_.resize(n_reg_, n_group_);
    flux_old_.resize(n_reg_, n_group_);
    vol_.resize(n_reg_);
    groups_ = Range(n_group_);

    // Set the mesh volumes. Same as the pin volumes
    int ipin = 0;
    for (auto &pin : mesh_) {
        int i   = mesh_.coarse_cell(mesh_.pin_position(ipin));
        vol_[i] = pin->vol();
        ipin++;
    }

    // Make sure we have input from the XML
    if (input.empty()) {
        throw EXCEPT("No input specified to initialize Sn sweeper.");
    }

    // Parse the number of inner iterations
    int int_in = input.attribute("n_inner").as_int(-1);
    if (int_in < 0) {
        throw EXCEPT("Invalid number of inner iterations specified "
                     "(n_inner).");
    }
    n_inner_ = int_in;

    // Try to read boundary update option
    if (!input.attribute("boundary_update").empty()) {
        std::string in_string = input.attribute("boundary_update").value();
        sanitize(in_string);

        if ((in_string == "gs") || (in_string == "gauss-seidel")) {
            gs_boundary_ = true;
        }
        else if ((in_string == "jacobi") || (in_string == "j")) {
            gs_boundary_ = false;
        }
        else {
            throw EXCEPT("Unrecognized option for BC update!");
        }
    }
    // For now, the BC doesnt support parallel boundary updates, so
    // disable Gauss-Seidel boundary update if we are using multiple
    // threads.
    /// \todo Add support for multi-threaded G-S boundary update in Sn
    LogScreen << omp_get_max_threads() << " " << gs_boundary_ << std::endl;
    if ((omp_get_max_threads() > 1) && gs_boundary_) {
        gs_boundary_ = false;
        Warn("Disabling Gauss-Seidel boundary update "
             "in parallel Sn");
    }

    timer_.toc();
    timer_init_.toc();

    return;
} // SnSweeper::SnSweeper( input, mesh )

/**
 * This just decides what method should be used to update the incoming flux
 * and instantiates an appropriate lambda function to carry out the update
 * for a single face. The lambda is then used to specialize the \ref
 * update_incoming_generic funcion template, which does most of the real
 * work.
 */
void SnSweeper::update_incoming_flux()
{
    // This should only be called in situations where coarse data should be
    // available
    assert(coarse_data_);

    // Short circuit if incoming flux update is desabled explicitly
    if (!do_incoming_update_) {
        return;
    }

    if (coarse_data_->has_old_partial()) {
        // There are old partial currents on the coarse data object, update
        // incomming flux based on the ratio of new to old
        auto update = [&](real_t in, int is, int g, int i) {
            real_t part = 2.0 * (coarse_data_->partial_current(is, g)[0] +
                                 coarse_data_->partial_current(is, g)[1]);
            real_t part_old =
                2.0 * (coarse_data_->partial_current_old(is, g)[0] +
                       coarse_data_->partial_current_old(is, g)[1]);

            if (part_old > 0.0) {
                real_t r = part / part_old;
                return in * r;
            }
            else {
                return in;
            }
        };
        update_incoming_generic<decltype(update)>(update);
    }
    else {
        // We dont have old partial currents to work with (probably because
        // this is the first iteration). Set the incomming flux to reflect
        // the partial current directly
        auto update = [&](real_t in, int is, int g, int i) {
            return (2.0 * RFPI * (coarse_data_->partial_current(is, g)[i] +
                                  coarse_data_->partial_current(is, g)[i]));
        };

        update_incoming_generic<decltype(update)>(update);
    }
}

/**
 * Watch out, this is potentially brittle, since it assumes parity between
 * the mesh regions and XS Mesh regions.
 */
ArrayB3 SnSweeper::pin_powers() const
{
    ArrayB3 powers(mesh_.nz(), mesh_.ny(), mesh_.nx());
    powers = 0.0;

    for (int ireg = 0; ireg < n_reg_; ireg++) {
        auto pos                = mesh_.coarse_position(ireg);
        const XSMeshRegion &xsr = (*xs_mesh_)[ireg];
        assert(xsr.reg().size() == 1);
        assert(xsr.reg()[0] == ireg);
        for (int ig = 0; ig < n_group_; ig++) {
            real_t p = vol_[ireg] * flux_(ireg, ig) * xsr.xsmacf(ig);
            powers(pos.z, pos.y, pos.x) += p;
        }
    }

    Normalize(powers.begin(), powers.end());

    return powers;
}

void SnSweeper::check_balance(int group) const
{
    if (!coarse_data_) {
        throw EXCEPT("No coarse data. Need it to look at currents.");
    }
    for (size_t icell = 0; icell < mesh_.n_pin(); icell++) {
        real_t b = 0.0;

        // Current
        b -= coarse_data_->current(mesh_.coarse_surf(icell, Surface::EAST),
                                   group) *
             mesh_.coarse_area(icell, Surface::EAST);
        b -= coarse_data_->current(mesh_.coarse_surf(icell, Surface::NORTH),
                                   group) *
             mesh_.coarse_area(icell, Surface::NORTH);
        b -= coarse_data_->current(mesh_.coarse_surf(icell, Surface::TOP),
                                   group) *
             mesh_.coarse_area(icell, Surface::TOP);
        b += coarse_data_->current(mesh_.coarse_surf(icell, Surface::WEST),
                                   group) *
             mesh_.coarse_area(icell, Surface::WEST);
        b += coarse_data_->current(mesh_.coarse_surf(icell, Surface::SOUTH),
                                   group) *
             mesh_.coarse_area(icell, Surface::SOUTH);
        b += coarse_data_->current(mesh_.coarse_surf(icell, Surface::BOTTOM),
                                   group) *
             mesh_.coarse_area(icell, Surface::BOTTOM);

        // Source
        b += (*source_)[icell] * vol_[icell];

        // Internal removal
        b -=
            flux_1g_(icell) * (*xs_mesh_)[icell].xsmacrm()[group] * vol_[icell];

        std::cout << "Cell balance: " << b << std::endl;
    }
    std::cout << std::endl;
}

void SnSweeper::output(H5Node &node) const
{
    auto dims = mesh_.dimensions();
    std::reverse(dims.begin(), dims.end());

    // Make a group in the file to store the flux
    node.create_group("flux");

    ArrayB2 flux = this->get_pin_flux();
    Normalize(flux.begin(), flux.end());

    for (int ig = 0; ig < n_group_; ig++) {
        std::stringstream setname;
        setname << "flux/" << std::setfill('0') << std::setw(3) << ig + 1;

        ArrayB1 flux_1g = flux(blitz::Range::all(), ig);

        node.write(setname.str(), flux_1g.begin(), flux_1g.end(), dims);
    }

    node.write("pin_powers", this->pin_powers());
    ang_quad_.output(node);

    LogFile << "Sn Sweeper:" << std::endl;

    LogFile << "Boundary update: ";
    if (gs_boundary_) {
        LogFile << "Gauss-Seidel" << std::endl;
    }
    else {
        LogFile << "Jacobi" << std::endl;
    }

    LogFile << std::endl;

    xs_mesh_->output(node);
    return;
}
} // namespace sn
} // namespace moc
