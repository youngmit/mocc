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

#include "moc_sweeper.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/string_utils.hpp"
#include "util/utils.hpp"
#include "util/validate_input.hpp"
#include "moc_current_worker.hpp"

namespace {
/**
 * \brief Return the appropriate sizing values for construcing a \ref
 * mocc::BoundaryCondition.
 *
 * This exists so that it may be used to construct the boundary condition
 * members of the \ref mocc::moc::MoCSweeper from the initializer list.
 */
std::vector<mocc::BC_Size_t> bc_size_helper(const mocc::moc::RayData &rays)
{
    std::vector<mocc::BC_Size_t> bc_dims(rays.ang_quad().ndir_oct() * 4);
    for (int iang = 0; iang < (int)rays.begin()->size(); iang++) {
        int iang1      = iang;
        int iang2      = rays.ang_quad().reverse(iang);
        int nx         = rays.ny(iang);
        int ny         = rays.nx(iang);
        bc_dims[iang1] = {{nx, ny, 0}};
        bc_dims[iang2] = {{nx, ny, 0}};
    }

    assert(bc_dims.size() == bc_dims.capacity());
    return bc_dims;
}

const std::vector<std::string> recognized_attributes = {
    "type",         "update_incoming", "n_inner",
    "dump_rays",    "boundary_update", "tl_splitting",
    "dump_fsr_flux"};
}

namespace mocc {
namespace moc {
MoCSweeper::MoCSweeper(const pugi::xml_node &input, const CoreMesh &mesh)
    : TransportSweeper(input, mesh, MeshTreatment::PLANE),
      timer_(RootTimer.new_timer("MoC Sweeper", true)),
      timer_init_(timer_.new_timer("Initialization", true)),
      timer_sweep_(timer_.new_timer("Sweep")),
      mesh_(mesh),
      rays_(input.child("rays"), ang_quad_, mesh),
      boundary_(mesh.nz(),
                BoundaryCondition(n_group_, ang_quad_, mesh_.boundary(),
                                  bc_size_helper(rays_))),
      boundary_out_(mesh.nz(), BoundaryCondition(1, ang_quad_, mesh_.boundary(),
                                                 bc_size_helper(rays_))),
      xstr_(xs_mesh_.get()),
      flux_1g_(),
      subplane_(mesh.subplane()),
      subplane_bounds_(),
      bc_type_(mesh_.boundary()),
      dump_rays_(false),
      dump_fsr_flux_(false),
      gauss_seidel_boundary_(true),
      allow_splitting_(false)
{
    LogFile << "Constructing a base MoC sweeper" << std::endl;

    validate_input(input, recognized_attributes);

    // Make sure we have input from the XML
    if (input.empty()) {
        throw EXCEPT("No input specified to initialize MoC sweeper.");
    }

    // Parse the number of inner iterations
    int int_in = input.attribute("n_inner").as_int(-1);
    if (int_in < 0) {
        throw EXCEPT("Invalid number of inner iterations specified "
                     "(n_inner).");
    }
    n_inner_ = int_in;

    // Parse the output options
    dump_rays_     = input.attribute("dump_rays").as_bool(false);
    dump_fsr_flux_ = input.attribute("dump_fsr_flux").as_bool(false);

    // Determine boundary update technique
    gauss_seidel_boundary_ = true;
    if (!input.attribute("boundary_update").empty()) {
        std::string in_string = input.attribute("boundary_update").value();
        sanitize(in_string);
        if (in_string == "jacobi" || in_string == "j") {
            gauss_seidel_boundary_ = false;
        } else if (in_string == "gs") {
        } else {
            throw EXCEPT("Unrecognized boundary update option.");
        }
    }

    // Parse TL source splitting setting
    allow_splitting_ = input.attribute("tl_splitting").as_bool(false);
    if (allow_splitting_) {
        split_.resize(n_reg_);
    }

    // Sanity-check the subplane parameters. We will operate on the assumption
    // for now that all planes in a macroplane are not only geometrically
    // identical, but completely so. For anyone interested in doing de-cusping,
    // this will need to change
    for (const auto &assembly : mesh_.core()) {
        int ip = 0;
        for (const auto &mac_size : subplane_) {
            int lat_id = (*assembly)[ip].id();
            for (int in_mac_plane = 1; in_mac_plane < mac_size;
                 in_mac_plane++) {
                if ((*assembly)[ip + in_mac_plane].id() != lat_id) {
                    throw EXCEPT(
                        "All lattices in a macroplane must be the same.");
                }
            }
            ip += mac_size;
        }
    }

    // Construct the entries in subplane_bounds_ these should be the the
    // accumulation of the number of planes in each MoC plane, such that taking
    // std::distance(..., std::lower_bound(..., iz)) for a given axial index iz
    // should return the appropriate MoC plane index
    subplane_bounds_.reserve(subplane_.size());
    int prev = 0;
    for (const auto np : subplane_) {
        subplane_bounds_.push_back(prev + np);
        prev += np;
    }

    // set up the array of unique plane ids and number of FSRs per moc plane. We
    // have already checked that all of the planes are entirely the same, so we
    // know they must be geometrically identical as well.
    macroplane_unique_ids_.reserve(subplane_.size());
    first_reg_macroplane_.reserve(subplane_.size());
    first_reg_macroplane_.push_back(0);
    nreg_plane_.reserve(subplane_.size());
    int iz = 0;
    for (const auto &sub_size : subplane_) {
        macroplane_unique_ids_.push_back(mesh_.unique_plane_ids()[iz]);
        nreg_plane_.push_back(
            mesh_.unique_plane(macroplane_unique_ids_.back()).n_reg());
        iz += sub_size;
    }

    for (const auto &plane_id : macroplane_unique_ids_) {
        first_reg_macroplane_.push_back(first_reg_macroplane_.back() +
                                        mesh_.unique_plane(plane_id).n_reg());
    }
    first_reg_macroplane_.pop_back();

    if (dump_rays_) {
        std::ofstream rayfile("rays.py");
        rayfile << rays_ << std::endl;
    }

    // Replace the angular quadrature with the modularized version
    ang_quad_ = rays_.ang_quad();

    timer_init_.toc();
    timer_.toc();

    return;
} // MoCSweeper( input, mesh )

void MoCSweeper::sweep(int group)
{
    assert(source_);

    timer_.tic();
    timer_sweep_.tic();

    // Expand the cross sections, and perform splitting if necessary
    xstr_.expand(group, split_);

    flux_1g_.reference(flux_(blitz::Range::all(), group));

    // Perform inner iterations
    for (unsigned int inner = 0; inner < n_inner_; inner++) {
        // update the self-scattering source
        if (inner ==0 && !(source_->get_has_external())) {
            source_->self_scatter_for_MMS(group, xstr_.xs());

            //print out the source and take a look at the result.
            std::cout << source_->n_reg() << std::endl;
            // std::cout << source_->get_source_1g_with_self_scat(0) << std::endl;

            std::string filename = "group_" + std::to_string(group+1) + "_source.txt";
            std::ofstream myfile;
            myfile.open (filename);
            myfile << source_->get_source_1g_with_self_scat(0);
            myfile.close();

            filename = "group_" + std::to_string(group+1) + "_xstr.txt";
            // std::ofstream myfile;
            myfile.open (filename);
            myfile << xstr_.xs();
            myfile.close();
            continue;
        }

        source_->self_scatter(group, xstr_.xs());

        // Perform the stock sweep unless we are on the last outer and have
        // a CoarseData object.
        if (inner == n_inner_ - 1 && coarse_data_) {
            // Wipe out the existing currents (only on X- and Y-normal
            // faces)
            coarse_data_->zero_data_radial(group);

            moc::Current cw(coarse_data_, &mesh_);
            this->sweep1g(group, cw);
            coarse_data_->set_has_radial_data(true);
        } else {
            moc::NoCurrent cw(coarse_data_, &mesh_);
            this->sweep1g(group, cw);
        }
    }

    timer_.toc();
    timer_sweep_.toc();
    return;
} // sweep( group )

/**
 * For now, this doesn't do anything remotely intelligent about the initial
 * guess for the scalar and angular flux values and just sets them to unity
 * and 1/4PI, respectively. At some point it might be useful to solve an IHM
 * problem and use at least the spectrum. In reality, it'd only cut down on
 * the initial CMFD iterations.
 */
void MoCSweeper::initialize()
{
    real_t val = 1.0; //0.5; //1.0;

    // Set the flux on the coarse mesh
    if (coarse_data_) {
        coarse_data_->flux = val;
    }
    // There are better ways to do this, but for now, just start with 1.0
    flux_     = val;
    flux_old_ = val;

    // Walk through the boundary conditions and initialize them the 1/4pi
    real_t bound_val = val / FPI;
    for (auto &boundary : boundary_) {
        boundary.initialize_scalar(bound_val);
    }

    return;
} // initialize()

void MoCSweeper::update_incoming_flux()
{
    assert(coarse_data_);

    // Short circuit if explicitly disabled
    if (!do_incoming_update_) {
        return;
    }

    if (coarse_data_->has_old_partial()) {
        auto update = [&](real_t in, int is, int ig) {
            real_t part = 2.0 * (coarse_data_->partial_current(is, ig)[0] +
                                 coarse_data_->partial_current(is, ig)[1]);
            real_t part_old =
                2.0 * (coarse_data_->partial_current_old(is, ig)[0] +
                       coarse_data_->partial_current_old(is, ig)[1]);

            if (part_old > 0.0) {
                real_t r = part / part_old;
                return in * r;
            } else {
                return in;
            }
        };
        update_incoming_generic<decltype(update)>(update);
    } else {
        auto update = [&](real_t in, int is, int ig) {
            real_t out =
                (2.0 * RFPI * (coarse_data_->partial_current(is, ig)[0] +
                               coarse_data_->partial_current(is, ig)[1]));
            return out;
        };
        update_incoming_generic<decltype(update)>(update);
    }

    return;
}

void MoCSweeper::get_pin_flux_1g(int group, ArrayB1 &flux,
                                 MeshTreatment treatment) const
{
    assert((int)flux.size() == mesh_.n_reg(treatment));
    /// \todo Put this back in when we address index ordering
    /// assert(flux.isStorageContiguous());
    flux = 0.0;

    switch (treatment) {
    case MeshTreatment::PIN_PLANE: {
        int ireg    = 0;
        int implane = 0;
        for (const auto &mplane : mesh_.macroplanes()) {
            int ipin = 0;
            for (const auto mpin : mplane) {
                real_t v        = 0.0;
                real_t pin_flux = 0.0;
                for (int ir = 0; ir < mpin->n_reg(); ir++) {
                    v += vol_[ireg];
                    pin_flux += flux_(ireg, group) * vol_[ireg];
                    ireg++;
                }
                pin_flux /= v;

                Position pos = mesh_.pin_position(ipin);
                pos.z        = implane;
                int i        = mesh_.coarse_cell(pos);
                flux(i) += pin_flux;

                ipin++;
            }
            implane++;
        }
    } break;

    case MeshTreatment::PIN: {
        int ireg = 0;
        for (const auto &mplane : mesh_.macroplanes()) {
            int ipin = 0;
            for (const auto mpin : mplane) {
                real_t v        = 0.0;
                real_t pin_flux = 0.0;
                for (int ir = 0; ir < mpin->n_reg(); ir++) {
                    v += vol_[ireg];
                    pin_flux += flux_(ireg, group) * vol_[ireg];
                    ireg++;
                }
                pin_flux /= v;

                Position pos = mesh_.pin_position(ipin);
                for (int iz = mplane.iz_min; iz <= mplane.iz_max; iz++) {
                    pos.z = iz;
                    int i = mesh_.coarse_cell(pos);
                    flux(i) += pin_flux;
                }

                ipin++;
            }
        }
    } break;
    default:
        throw EXCEPT("Unsupported mesh treatment requested");
    }

    return;
}

real_t MoCSweeper::set_pin_flux_1g(int group, const ArrayB1 &pin_flux,
                                   MeshTreatment treatment)
{
    assert((int)pin_flux.size() == mesh_.n_reg(treatment));

    ArrayB1 plane_pin_flux(mesh_.nx() * mesh_.ny() * subplane_.size());

    // Check for setting any of the pin fluxes to zero. This can cause lots of
    // issues down the line.
    for(auto v: pin_flux) {
        if(v <= 0.0 ) {
            std::stringstream msg;
            msg << "Negative or zero input flux: " << v;
            throw EXCEPT(msg.str());
        }
    }

    real_t resid = 0.0;

    // Watch out: we use fall-through here, since once we have a
    // macroplane-homogenized flux, we can use the same logic.
    switch (treatment) {
    case MeshTreatment::PIN:
        plane_pin_flux = 0.0;
        // homogenize the passed-in pin flux to the coarser axial mesh.
        for (int i = 0; i < mesh_.n_pin(); i++) {
            Position pos = mesh_.coarse_position(i);
            int iz       = pos.z;
            pos.z        = this->moc_plane_index(pos.z);
            plane_pin_flux(mesh_.coarse_cell(pos)) +=
                pin_flux(i) * mesh_.dz(iz);
        }
        for (unsigned i = 0; i < plane_pin_flux.size(); i++) {
            int iz = i / (mesh_.nx() * mesh_.ny());
            plane_pin_flux(i) /= mesh_.macroplane_heights()[iz];
        }
    case MeshTreatment::PIN_PLANE: {
        plane_pin_flux = pin_flux;
        int iz       = 0;
        int ireg     = 0;
        for (const auto &mplane : mesh_.macroplanes()) {
            int ipin = 0;
            for (const auto &pin : mplane) {
                Position pos   = mesh_.pin_position(ipin);
                pos.z          = iz;
                int i_coarse   = mesh_.coarse_cell(pos);
                real_t fm_flux = 0.0;
                for (const auto area : pin->areas()) {
                    fm_flux += flux_(ireg, group) * area;
                    ireg++;
                }
                fm_flux /= pin->area();
                real_t e = plane_pin_flux(i_coarse) - fm_flux;
                real_t f = plane_pin_flux(i_coarse) / fm_flux;
                ireg -= pin->n_reg();

                for (int ir = 0; ir < pin->n_reg(); ir++) {
                    flux_(ireg, group) *= f;
                    ireg++;
                }

                resid += e * e;
                ipin++;
            }
            iz++;
        }
    } break;
    default:
        // i suppose there is no reason we couldn't use TRUE...
        throw EXCEPT("Unsupported mesh treatment used");
    }

    return std::sqrt(resid) / plane_pin_flux.size();
} // set_pin_flux_1f( group, pin_flux )

void MoCSweeper::apply_transverse_leakage(int group, const ArrayB1 &tl)
{
    assert((int)tl.size() == n_reg_);

    flux_1g_.reference(flux_(blitz::Range::all(), group));

    /// \todo for now, this is using a pretty invasive direct access the the
    /// source. Might be good to do as a call to auxiliary() instead
    if (allow_splitting_) {
        int n_split = 0;
        split_      = 0.0;
        for (int ireg = 0; ireg < n_reg_; ireg++) {
            real_t s = (*source_)[ireg] + tl(ireg);
            if (s < 0.0) {
                if (flux_1g_(ireg) < 0.0) {
                    std::cout << ireg << " " << flux_1g_(ireg) << std::endl;
                    throw EXCEPT("Negative flux when splitting");
                }
                n_split++;
                split_(ireg)     = -s / flux_1g_(ireg);
                (*source_)[ireg] = 0.0;
            } else {
                (*source_)[ireg] = s;
            }
        }

        if (n_split > 0) {
            LogFile << "Split " << n_split << " region sources" << std::endl;
        }
    } else {
        for (int ireg = 0; ireg < n_reg_; ireg++) {
            real_t s         = (*source_)[ireg] + tl(ireg);
            (*source_)[ireg] = s;
        }
    }

    return;
} // apply_transverse_leakage( tl )

/**
 * \brief Check for the balance of neutrons within each pin cell.
 *
 * \todo Make sure this is valid in the presence of source splitting
 */
void MoCSweeper::check_balance(int group) const
{
    ArrayB1 b(mesh_.n_pin());
    b = 0.0;

    // Get the removal cross section in a nice format
    ArrayB1 xsrm(n_reg_);
    xsrm = 0.0;
    for (auto &xsr : *xs_mesh_) {
        real_t rm = xsr.xsmacrm(group);
        for (auto &ireg : xsr.reg()) {
            xsrm(ireg) = rm;
        }
    }

    ArrayB1 current_1g = coarse_data_->current(blitz::Range::all(), group);

    int ipin = 0;
    int ireg = 0;
    for (const auto pin : mesh_) {
        int icell = mesh_.coarse_cell(mesh_.pin_position(ipin));
        real_t bi = 0.0;

        for (int ireg_pin = 0; ireg_pin < pin->n_reg(); ireg_pin++) {
            bi -= flux_(ireg, group) * vol_[ireg] * xsrm(ireg);
            bi += (*source_)[ireg] * vol_[ireg];
            ireg++;
        }

        // Current
        bi -= current_1g(mesh_.coarse_surf(icell, Surface::EAST)) *
              mesh_.coarse_area(icell, Surface::EAST);
        bi -= current_1g(mesh_.coarse_surf(icell, Surface::NORTH)) *
              mesh_.coarse_area(icell, Surface::NORTH);
        bi -= current_1g(mesh_.coarse_surf(icell, Surface::TOP)) *
              mesh_.coarse_area(icell, Surface::TOP);
        bi += current_1g(mesh_.coarse_surf(icell, Surface::WEST)) *
              mesh_.coarse_area(icell, Surface::WEST);
        bi += current_1g(mesh_.coarse_surf(icell, Surface::SOUTH)) *
              mesh_.coarse_area(icell, Surface::SOUTH);
        bi += current_1g(mesh_.coarse_surf(icell, Surface::BOTTOM)) *
              mesh_.coarse_area(icell, Surface::BOTTOM);

        b(icell) = bi;
        ipin++;
    }

    std::cout << "MoC cell balance:" << std::endl;
    for (auto v : b) {
        std::cout << v << std::endl;
    }

    return;
}

void MoCSweeper::output(H5Node &node) const
{
    // Get core dimensions from the mesh
    VecI dims = mesh_.dimensions();
    std::reverse(dims.begin(), dims.end());

    // Make a group in the file to store the flux
    node.create_group("flux");

    ArrayB2 flux = this->get_pin_flux();
    Normalize(flux.begin(), flux.end());

    // Make a group in the file to store the fsr_flux if dump_fsr_flux is trur
    if (dump_fsr_flux_) {
        VecI dims = {1, (int)(mesh_.n_reg(MeshTreatment::PLANE))};
        std::reverse(dims.begin(), dims.end());

        node.create_group("fsr_flux");

        for (int ig = 0; ig < n_group_; ig++) {
            std::stringstream setname;
            setname << "fsr_flux/" << std::setfill('0') << std::setw(3)
                    << ig + 1;

            ArrayB1 flux_1g = flux_(blitz::Range::all(), ig);

            node.write(setname.str(), flux_1g.begin(), flux_1g.end(), dims);
        }
    }

    LogFile << "Boundary update: ";
    if (gauss_seidel_boundary_) {
        LogFile << "Gauss-Seidel" << std::endl;
    } else {
        LogFile << "Jacobi" << std::endl;
    }

    for (int ig = 0; ig < n_group_; ig++) {
        std::stringstream setname;
        setname << "flux/" << std::setfill('0') << std::setw(3) << ig + 1;

        ArrayB1 flux_1g = flux(blitz::Range::all(), ig);
        node.write(setname.str(), flux_1g.begin(), flux_1g.end(), dims);
    }

    // Pin powers
    node.write("pin_powers", this->pin_powers());

    ang_quad_.output(node);

    return;
}
}
}
