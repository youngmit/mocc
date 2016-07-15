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

#include "core_mesh.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/string_utils.hpp"

namespace mocc {
CoreMesh::CoreMesh(const pugi::xml_node &input)
    : pin_meshes_(ParsePinMeshes(input)),
      mat_lib_(input.child("material_lib")),
      pins_(ParsePins(input, pin_meshes_, mat_lib_)),
      lattices_(ParseLattices(input, pins_)),
      assemblies_(ParseAssemblies(input, lattices_)),
      core_(ParseCore(input, assemblies_)),
      subplane_(core_.front().subplane())
{
    LogScreen << "Building core mesh... " << std::endl;

    nx_           = core_.npin_x();
    ny_           = core_.npin_y();
    nz_           = core_.nz();
    n_surf_plane_ = (nx_ + 1) * ny_ + (ny_ + 1) * nx_ + nx_ * ny_;
    nasy_         = core_.nasy();
    bc_           = core_.boundary();

    // Calculate the total core dimensions
    hx_ = 0.0;
    for (int ix = 0; ix < core_.nx(); ix++) {
        hx_ += core_.at(ix, 0).hx();
    }
    hy_ = 0.0;
    for (int iy = 0; iy < core_.ny(); iy++) {
        hy_ += core_.at(0, iy).hy();
    }

    // Determine the set of geometrically-unique axial planes and plane region
    // offsets
    std::vector<const Lattice *> plane_lattices;
    plane_lattices.reserve(core_.nx() * core_.ny());
    int plane_reg = 0;
    n_fuel_2d_    = 0;
    for (int iz = 0; iz < nz_; iz++) {
        first_reg_plane_.push_back(plane_reg);
        // Construct a temporary Plane object to test against the list of known
        // unique planes, or to insert if its new
        for (const auto &assembly : core_) {
            const Lattice &lattice = (*assembly)[iz];
            plane_lattices.push_back(&lattice);
            plane_reg += (*assembly)[iz].n_reg();
            // add the contained pins to the overall list of pins in the
            // CoreMesh
            for (const auto &pin : lattice) {
                core_pins_.push_back(pin);
            }
        }
        Plane temp_plane(plane_lattices, core_.nx(), core_.ny());
        n_fuel_2d_ = std::max(n_fuel_2d_, temp_plane.n_fuel());

        // Check against current list of unique planes
        auto match_plane =
            std::find_if(planes_.begin(), planes_.end(), [&](const Plane &p) {
                return p.geometrically_equivalent(temp_plane);
            });
        int unique_id = std::distance(planes_.begin(), match_plane);
        if (match_plane == planes_.end()) {
            // This plane is thus far unique. Add it to planes_
            planes_.push_back(temp_plane);
            first_unique_.push_back(iz);
        } else {
            // This plane is geometrically identical to another plane we have
            // already visited. do nothing
        }
        unique_plane_ids_.push_back(unique_id);
        plane_lattices.clear();
    } // Unique plane search
    LogFile << "Unique plane search done" << std::endl;

    // Put together the list of pin boundaries. For now we are treating them
    // as independent of axial plane
    x_vec_.push_back(0.0);
    real_t h_prev = 0.0;
    for (int ilatx = 0; ilatx < core_.nx(); ilatx++) {
        const Assembly *asy = &core_.at(ilatx, 0);
        const Lattice *lat  = &((*asy)[0]);
        for (auto h : lat->hx_vec()) {
            dx_vec_.push_back(h);
            x_vec_.push_back(h + h_prev);
            h_prev += h;
        }
    }
    y_vec_.push_back(0.0);
    h_prev = 0.0;
    for (int ilaty = 0; ilaty < core_.ny(); ilaty++) {
        const Assembly *asy = &core_.at(0, ilaty);
        const Lattice *lat  = &((*asy)[0]);
        for (auto h : lat->hy_vec()) {
            dy_vec_.push_back(h);
            y_vec_.push_back(h + h_prev);
            h_prev += h;
        }
    }

    // Form lines for internal pin boundaries and domain bounding box
    for (auto xi = x_vec_.cbegin() + 1; xi != x_vec_.cend() - 1; ++xi) {
        lines_.push_back(Line(Point2(*xi, 0.0), Point2(*xi, hx_)));
    }
    for (auto yi = y_vec_.cbegin() + 1; yi != y_vec_.cend() - 1; ++yi) {
        lines_.push_back(Line(Point2(0.0, *yi), Point2(hx_, *yi)));
    }
    bounding_box_ = Box(Point2(0.0, 0.0), Point2(hx_, hy_));

    dz_vec_ = core_.dz();
    z_vec_.reserve(dz_vec_.size() + 1);
    hz_ = 0.0;
    z_vec_.push_back(hz_);
    for (const auto &dz : dz_vec_) {
        hz_ += dz;
        z_vec_.push_back(hz_);
    }

    // Coarse mesh volumes.
    /**
     * \todo Really need to clean up \ref Mesh and \ref CoreMesh
     * construction. The data members of the base \ref Mesh class should
     * really be initialized by the Mesh class itself. Its hard to know the
     * necessary data to pass to the Mesh constructor before parsing the
     * core first. Going to have to come up with something.
     */
    coarse_vol_ = VecF(this->n_pin());
    for (int i = 0; i < this->n_pin(); i++) {
        auto pos       = this->coarse_position(i);
        coarse_vol_[i] = dx_vec_[pos.x] * dy_vec_[pos.y] * dz_vec_[pos.z];
    }

    // Add up the number of regions and XS regions in the entire problem
    // geometry
    n_reg_   = 0;
    n_xsreg_ = 0;
    for (auto &a : core_.assemblies()) {
        n_reg_ += a->n_reg();
        n_xsreg_ += a->n_xsreg();
    }

    // Figure out the height of the macroplanes
    macroplane_heights_ = VecF(subplane_.size(), 0.0);
    int iz              = 0;
    int isub            = 0;
    for (const auto &np : subplane_) {
        for (int i = 0; i < np; i++) {
            macroplane_heights_[isub] += dz_vec_[iz];
            iz++;
        }
        isub++;
    }

    // Figure out the macroplanes in general
    macroplanes_.reserve(subplane_.size());
    iz         = 0;
    int iplane = 0;
    for (const auto np : subplane_) {
        MacroPlane::PinIter_t stt = core_pins_.begin() + iz * nx_ * ny_;
        MacroPlane::PinIter_t stp = stt + nx_ * ny_;
        macroplanes_.emplace_back(&(planes_[unique_plane_ids_[iz]]), iz,
                                  iz + np - 1, macroplane_heights_[iplane], stt,
                                  stp);
        iz += np;
        iplane++;
    }

    // Set the macroplane indexes on the base Mesh type
    macroplane_index_.resize(nz_);
    iplane = 0;
    for (const auto &mplane : macroplanes_) {
        for (int iz = mplane.iz_min; iz <= mplane.iz_max; iz++) {
            macroplane_index_[iz] = iplane;
        }
        iplane++;
    }

    // calculate surface indices
    this->prepare_surfaces();

    LogScreen << "Done building Core Mesh." << std::endl;

    return;
} // constructor

CoreMesh::~CoreMesh()
{
    return;
}

const PinMeshTuple CoreMesh::get_pinmesh(Point2 &p, size_t iplane,
                                         int &first_reg) const
{
    assert((iplane >= 0) & (iplane < planes_.size()));

    // Locate the Position of the pin
    unsigned int ix = std::lower_bound(x_vec_.begin(), x_vec_.end(), p.x) -
                      x_vec_.begin() - 1;
    unsigned int iy = std::lower_bound(y_vec_.begin(), y_vec_.end(), p.y) -
                      y_vec_.begin() - 1;

    Position pos(ix, iy, 0);

    return PinMeshTuple(pos, planes_[iplane].get_pinmesh(p, first_reg));
}

Position CoreMesh::pin_position(size_t ipin) const
{
    Position pos = planes_[0].pin_position(ipin % (nx_ * ny_));
    pos.z        = ipin / (nx_ * ny_);
    return pos;
}

Point2 CoreMesh::pin_origin(size_t ipin) const
{
    auto pos = this->pin_position(ipin);
    Point2 p;
    p.x = (x_vec_[pos.x] + x_vec_[pos.x + 1]) * 0.5;
    p.y = (y_vec_[pos.y] + y_vec_[pos.y + 1]) * 0.5;

    return p;
}

CoreMesh::LocationInfo CoreMesh::get_location_info(Point3 p,
                                                   Direction dir) const
{
    LocationInfo info;

    // Locate the pin Position
    int ix = std::distance(x_vec_.begin(),
                           std::lower_bound(x_vec_.begin(), x_vec_.end(), p.x));
    info.pos.x = ix - 1;

    int iy = std::distance(y_vec_.begin(),
                           std::lower_bound(y_vec_.begin(), y_vec_.end(), p.y));
    info.pos.y = iy - 1;

    info.pos.z = this->plane_index(p.z, dir.oz);

    // Make a 2D copy of the 3D point to convert to pin-local
    info.local_point = p.to_2d();

    info.reg_offset = first_reg_plane(info.pos.z);

    Point2 pin_origin = p.to_2d();
    int plane_index   = unique_plane_ids_[info.pos.z];
    info.pm =
        planes_[plane_index].get_pinmesh(pin_origin, info.reg_offset, dir);
    info.local_point -= pin_origin;

    info.pin_boundary = {
        {Point3(x_vec_[info.pos.x], y_vec_[info.pos.y], z_vec_[info.pos.z]),
         Point3(x_vec_[info.pos.x + 1], y_vec_[info.pos.y + 1],
                z_vec_[info.pos.z + 1])}};

    return info;
}

std::ostream &operator<<(std::ostream &os, const CoreMesh &mesh)
{
    os << "Boundary conditions: " << std::endl;
    for (int ib = 0; ib < 6; ib++) {
        os << (Surface)ib << ":\t" << mesh.bc_[ib] << std::endl;
    }
    os << std::endl;

    os << "Mesh X Pitches:" << std::endl;
    for (auto v : mesh.dx_vec_) {
        os << v << std::endl;
    }
    os << std::endl;

    os << "Mesh Y Pitches:" << std::endl;
    for (auto v : mesh.dy_vec_) {
        os << v << std::endl;
    }
    os << std::endl;

    os << "Mesh Z Pitches:" << std::endl;
    for (auto v : mesh.dz_vec_) {
        os << v << std::endl;
    }
    os << std::endl;

    os << "Pin Meshes: " << std::endl;
    for (const auto &pm : mesh.pin_meshes_) {
        os << "Mesh ID: " << pm.first << std::endl;
        os << *(pm.second) << std::endl;
        os << std::endl;
    }

    return os;
}
}
