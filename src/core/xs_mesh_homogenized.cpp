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

#include "xs_mesh_homogenized.hpp"

#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include "pugixml.hpp"
#include "util/files.hpp"
#include "util/h5file.hpp"
#include "util/string_utils.hpp"

namespace mocc {
XSMeshHomogenized::XSMeshHomogenized(const CoreMesh &mesh)
    : XSMesh(mesh, MeshTreatment::PIN), mesh_(mesh), flux_(nullptr)
{
    // Set up the non-xs part of the xs mesh
    eubounds_ = mesh_.mat_lib().g_bounds();
    ng_       = eubounds_.size();

    int n_xsreg = std::accumulate(
        mesh_.macroplanes().begin(), mesh_.macroplanes().end(), 0,
        [](int n, const MacroPlane p) { return n + p.plane->n_pin(); });

    // Allocate space to store the cross sections
    this->allocate_xs(n_xsreg, ng_);

    // Set up the regions
    //
    // It is assumed that this type of XS mesh will be used to provide
    // cross sections to a MeshTreatment::PIN type sweeper, possibly using
    // flux-weighted homogenization with flux from a MeshTreatment::PLANE-type
    // sweeper. The mapping from XS mesh region to sweeper region will therefore
    // be from a homogenized pin at the macroplane level to all of the regions
    // in a MeshTreatment::PIN-type mesh within the corresponding macroplane.
    regions_.reserve(n_xsreg);
    n_reg_expanded_ = 0;
    int ixsreg = 0;
    for (const auto &mplane : mesh_.macroplanes()) {
        VecI ireg;
        ireg.reserve(mplane.iz_max - mplane.iz_min + 1);
        for (unsigned ipin = 0; ipin < mplane.size(); ++ipin) {
            Position pos = mesh_.pin_position(ipin);
            for (int iz = mplane.iz_min; iz <= mplane.iz_max; ++iz) {
                pos.z = iz;
                ireg.push_back(mesh_.coarse_cell(pos));
            }
            regions_.emplace_back(ireg, &xstr_(ixsreg, 0), &xsnf_(ixsreg, 0),
                                  &xsch_(ixsreg, 0), &xsf_(ixsreg, 0),
                                  &xsrm_(ixsreg, 0), ScatteringMatrix());

            n_reg_expanded_ += ireg.size();
            ireg.clear();
            ixsreg++;
        }
    }

    // Homogenize initial cross sections
    ixsreg = 0;
    for (const auto &mplane : mesh_.macroplanes()) {
        for (const auto &mpin : mplane) {
            this->homogenize_region(ixsreg, *mpin, regions_[ixsreg]);
            ixsreg++;
        }
    }
}

/**
 * This creates a "homogenized" cross-section mesh using data prescribed by
 * one or more HDF5 files. The allowance of multiple files is so that an
 * \ref XSMeshHomogenized object can be composed from various sources. For
 * example, imagine a 3-D mesh, where the cross sections are obtained from
 * several different 2-D solutions, and the cross sections for each plane
 * may be drawn from any of the 2-D results.
 *
 * Since the primary use-case for this functionality is for using planar,
 * fine-mesh MoC solutions to provide cross sections to a coarse-mesh Sn
 * solver, the data is specified along with the upper bound plane index to
 * which the cross sections are to be applied. The data itself is assumed to
 * be coming from a call to \ref XSMesh::output().
 */
XSMeshHomogenized::XSMeshHomogenized(const CoreMesh &mesh,
                                     const pugi::xml_node &input)
    : XSMeshHomogenized(mesh)
{
    if (input.child("data").empty()) {
        throw EXCEPT("No data found in input tag.");
    }

    int nreg_plane = mesh.nx() * mesh.ny();

    // First, validate the data tags. Make sure that they are the right size and
    // cover valid planes in the mesh.
    {
        int np = mesh_.macroplanes().size();
        std::vector<bool> plane_data(np, false);
        for (auto data = input.child("data"); data;
             data      = data.next_sibling("data")) {
            int top_plane = 0;
            int bot_plane = 0;
            if (!data.attribute("top_plane").empty()) {
                top_plane = data.attribute("top_plane").as_int();
            }
            if (!data.attribute("bottom_plane").empty()) {
                bot_plane = data.attribute("bottom_plane").as_int();
            }

            if ((bot_plane < 0) || (bot_plane >= np)) {
                throw EXCEPT("Invalid bottom_plane");
            }
            if ((top_plane < 0) || (top_plane >= np)) {
                throw EXCEPT("Invalid top_plane");
            }

            if (data.attribute("file").empty()) {
                throw EXCEPT("No file specified.");
            }

            // Check the dimensions of the contained XSMesh
            H5Node h5d(data.attribute("file").value(), H5Access::READ);
            auto dims = h5d.dimensions("/xsmesh/xstr/0");
            if (dims[2] != mesh_.nx()) {
                std::stringstream msg;
                msg << "Incorrect XS Mesh dimensions: " << dims[2] << " "
                    << mesh_.nx();
                throw EXCEPT(msg.str());
            }
            if (dims[1] != mesh_.ny()) {
                std::stringstream msg;
                msg << "Incorrect XS Mesh dimensions: " << dims[1] << " "
                    << mesh_.ny();
                throw EXCEPT("Incorrect XS Mesh dimensions");
            }
            if (dims[0] != 1) {
                throw EXCEPT("Data should only have one plane");
            }

            for (int ip = bot_plane; ip <= top_plane; ip++) {
                if (plane_data[ip]) {
                    std::stringstream msg;
                    msg << "Plane data is over-specified. Look at plane " << ip;
                    throw EXCEPT(msg.str());
                }
                plane_data[ip] = true;
            }
        }

        LogFile << "Data is being specified for the following macroplanes:"
                << std::endl;
        LogFile << print_range(plane_data) << std::endl;
    }

    // We are going to want a nice contiguous buffer for storing xs as they
    // come in from the file
    ArrayB1 tr_buf(nreg_plane);
    ArrayB1 nf_buf(nreg_plane);
    ArrayB1 ch_buf(nreg_plane);
    ArrayB1 kf_buf(nreg_plane);

    for (auto data = input.child("data"); data;
         data      = data.next_sibling("data")) {
        int top_plane = 0;
        int bot_plane = 0;
        if (!data.attribute("top_plane").empty()) {
            top_plane = data.attribute("top_plane").as_int();
        }
        if (!data.attribute("bottom_plane").empty()) {
            bot_plane = data.attribute("bottom_plane").as_int();
        }

        H5Node h5d(data.attribute("file").value(), H5Access::READ);

        h5d.read("//xsmesh/eubounds", eubounds_);

        // Get all the group data out to memory first
        for (int ig = 0; ig < (int)ng_; ig++) {
            {
                std::stringstream path;
                path << "/xsmesh/xstr/" << ig;
                h5d.read(path.str(), tr_buf);
            }
            {
                std::stringstream path;
                path << "/xsmesh/xsnf/" << ig;
                h5d.read(path.str(), nf_buf);
            }
            {
                std::stringstream path;
                path << "/xsmesh/xsch/" << ig;
                h5d.read(path.str(), ch_buf);
            }
            {
                std::stringstream path;
                path << "/xsmesh/xsf/" << ig;
                h5d.read(path.str(), kf_buf);
            }

            // Apply the above to the appropriate planes of the actual xs mesh
            // Since the data are output to the HDF 5 for plotting, they need to
            // be rearranged to the pin order
            for (int ip = bot_plane; ip <= top_plane; ip++) {
                const MacroPlane &mplane = mesh_.macroplanes()[ip];
                for (unsigned ipin = 0; ipin < mplane.size(); ipin++) {
                    int ixsreg   = nreg_plane * ip + ipin;
                    Position pos = mesh_.pin_position(ipin);
                    pos.z        = 0;
                    int icell    = mesh_.coarse_cell(pos);
                    xstr_(ixsreg, ig) = tr_buf(icell);
                    xsnf_(ixsreg, ig) = nf_buf(icell);
                    xsch_(ixsreg, ig) = ch_buf(icell);
                    xsf_(ixsreg, ig)  = kf_buf(icell);
                }
            }
        }

        // We dont try to plot the scattering cross sections in the same way
        // as we do the others, so this can be read in more naturally.
        ArrayB3 scat;
        h5d.read("/xsmesh/xssc", scat);

        // Storing the results of the read above into the appropriate arrays
        // should be sufficient for all cross sections except the scattering
        // matrices. Go through and store the scattering matrices on the xsmesh
        // regions too
        for (int ipin = 0; ipin < nreg_plane; ipin++) {
            ScatteringMatrix xssc(
                scat(ipin, blitz::Range::all(), blitz::Range::all()));
            for (int ip = bot_plane; ip <= top_plane; ip++) {
                int ixsreg                = ipin + nreg_plane * ip;
                regions_[ixsreg].xsmacsc_ = xssc;
            }
        }
    } // <data> tag loop

    return;
} // HDF5 constructor

/**
 * Update the XS mesh, incorporating a new estimate of the scalar flux.
 */
void XSMeshHomogenized::update()
{
    if (!flux_) {
        // no update needed if doing volume-weighted cross sections
        return;
    } else {
        // For now assume that the flux is coming from a PLANE-type sweeper
        assert(flux_->extent(0) == (int)mesh_.n_reg(MeshTreatment::PLANE));
    }
    int first_reg = 0;
    int ixsreg    = 0;
    for (const auto &mplane : mesh_.macroplanes()) {
        for (const auto &pin : mplane) {
            this->homogenize_region_flux(ixsreg, first_reg, *pin,
                                         regions_[ixsreg]);
            first_reg += pin->n_reg();
            ixsreg++;
        }
    }
    state_++;
    return;
}

void XSMeshHomogenized::homogenize_region(int i, const Pin &pin,
                                          XSMeshRegion &xsr) const
{
    VecI fsrs(1, i);
    VecF xstr(ng_, 0.0);
    VecF xsnf(ng_, 0.0);
    VecF xsf(ng_, 0.0);
    VecF xsch(ng_, 0.0);

    std::vector<VecF> scat(ng_, VecF(ng_, 0.0));

    auto mat_lib   = mesh_.mat_lib();
    auto &pin_mesh = pin.mesh();
    auto areas     = pin_mesh.areas();

    for (size_t ig = 0; ig < ng_; ig++) {
        int ireg     = 0;
        int ixsreg   = 0;
        real_t farea = 0.0;
        for (auto &mat_id : pin.mat_ids()) {
            auto mat                      = mat_lib.get_material_by_id(mat_id);
            const ScatteringRow &scat_row = mat.xssc().to(ig);
            int gmin                      = scat_row.min_g;
            int gmax                      = scat_row.max_g;
            real_t fsrc                   = 0.0;
            for (size_t igg = 0; igg < ng_; igg++) {
                fsrc += mat.xsnf(igg);
            }
            for (size_t i = 0; i < pin_mesh.n_fsrs(ixsreg); i++) {
                farea += areas[ireg] * fsrc;
                xstr[ig] += areas[ireg] * mat.xstr(ig);
                xsnf[ig] += areas[ireg] * mat.xsnf(ig);
                xsf[ig] += areas[ireg] * mat.xsf(ig);
                xsch[ig] += areas[ireg] * fsrc * mat.xsch(ig);

                for (int igg = gmin; igg <= gmax; igg++) {
                    scat[ig][igg] += scat_row.from[igg - gmin] * areas[ireg];
                }
                ireg++;
            }
            ixsreg++;
        }

        xstr[ig] /= pin.area();
        xsnf[ig] /= pin.area();
        xsf[ig] /= pin.area();
        if (farea > 0.0) {
            xsch[ig] /= farea;
        }

        for (auto &s : scat[ig]) {
            s /= pin.area();
        }
    }

    ScatteringMatrix scat_mat(scat);

    xsr.update(xstr, xsnf, xsch, xsf, scat_mat);

    return;
}

void XSMeshHomogenized::homogenize_region_flux(int i, int first_reg,
                                               const Pin &pin,
                                               XSMeshRegion &xsr) const
{
    assert(flux_);
    // Extract a reference to the flux array.
    const ArrayB2 flux = *flux_;

    // Set the FSRs to be one element, representing the coarse mesh index to
    // which this \ref XSMeshRegion belongs.
    VecI fsrs(1, i);
    VecF xstr(ng_, 0.0);
    VecF xsnf(ng_, 0.0);
    VecF xsf(ng_, 0.0);
    VecF xsch(ng_, 0.0);

    std::vector<VecF> scat(ng_, VecF(ng_, 0.0));

    auto mat_lib   = mesh_.mat_lib();
    auto &pin_mesh = pin.mesh();
    auto areas     = pin_mesh.areas();

    // Precompute the fission source in each region, since it is the
    // wieghting factor for chi
    VecF fs(pin_mesh.n_reg(), 0.0);
    {
        int ixsreg = 0;
        for (auto &mat_id : pin.mat_ids()) {
            auto mat = mat_lib.get_material_by_id(mat_id);
            for (size_t ig = 0; ig < ng_; ig++) {
                int ireg       = first_reg; // pin-local region index
                int ireg_local = 0;
                for (size_t i = 0; i < pin_mesh.n_fsrs(ixsreg); i++) {
                    fs[ireg_local] +=
                        mat.xsnf(ig) * flux(ireg, (int)ig) * areas[ireg_local];
                    ireg++;
                    ireg_local++;
                }
            }
            ixsreg++;
        }
    }

    real_t fs_sum = 0.0;
    for (auto &v : fs) {
        fs_sum += v;
    }

    for (size_t ig = 0; ig < ng_; ig++) {
        real_t fluxvolsum = 0.0;
        VecF scatsum(ng_, 0.0);
        int ireg       = first_reg; // global region index
        int ireg_local = 0;         // pin-local refion index
        int ixsreg     = 0;
        for (auto &mat_id : pin.mat_ids()) {
            auto mat                      = mat_lib.get_material_by_id(mat_id);
            const ScatteringRow &scat_row = mat.xssc().to(ig);
            size_t gmin                   = scat_row.min_g;
            size_t gmax                   = scat_row.max_g;
            for (size_t i = 0; i < pin_mesh.n_fsrs(ixsreg); i++) {
                real_t v      = areas[ireg_local];
                real_t flux_i = flux(ireg, (int)ig);
                fluxvolsum += v * flux_i;
                xstr[ig] += v * flux_i * mat.xstr(ig);
                xsnf[ig] += v * flux_i * mat.xsnf(ig);
                xsf[ig] += v * flux_i * mat.xsf(ig);
                xsch[ig] += fs[ireg_local] * mat.xsch(ig);

                for (size_t igg = 0; igg < ng_; igg++) {
                    real_t fluxgg = flux(ireg, (int)igg);
                    scatsum[igg] += fluxgg * v;
                    if ((igg >= gmin) && (igg <= gmax)) {
                        real_t scgg = scat_row.from[igg - gmin];
                        scat[ig][igg] += scgg * v * fluxgg;
                    }
                }
                ireg++;
                ireg_local++;
            }
            ixsreg++;
        }

        for (size_t igg = 0; igg < ng_; igg++) {
            if (scat[ig][igg] > 0.0) {
                scat[ig][igg] /= scatsum[igg];
            }
        }

        xstr[ig] /= fluxvolsum;
        xsnf[ig] /= fluxvolsum;
        xsf[ig] /= fluxvolsum;
        if (fs_sum > 0.0) {
            xsch[ig] /= fs_sum;
        }
    }

    ScatteringMatrix scat_mat(scat);

    xsr.update(xstr, xsnf, xsch, xsf, scat_mat);

    return;
}

void XSMeshHomogenized::output(H5Node &file) const
{
    auto xsm_g = file.create_group("xsmesh");
    xsm_g.create_group("xstr");
    xsm_g.create_group("xsnf");
    xsm_g.create_group("xsf");
    xsm_g.create_group("xsch");

    xsm_g.write("eubounds", eubounds_, VecI(1, ng_));

    auto d = mesh_.dimensions();
    std::reverse(d.begin(), d.end());
    d.front()     = mesh_.macroplanes().size();
    int per_plane = d[1] * d[2];

    for (size_t ig = 0; ig < ng_; ig++) {
        // Transport cross section
        VecF xstr(this->size(), 0.0);
        VecF xsnf(this->size(), 0.0);
        VecF xsf(this->size(), 0.0);
        VecF xsch(this->size(), 0.0);
        int i = 0;
        for (auto xsr : regions_) {
            Position pos = mesh_.pin_position(i);
            pos.z        = i / per_plane;
            int icell    = mesh_.coarse_cell(pos);
            xstr[icell]  = xsr.xsmactr(ig);
            xsnf[icell]  = xsr.xsmacnf(ig);
            xsf[icell]   = xsr.xsmacf(ig);
            xsch[icell]  = xsr.xsmacch(ig);
            i++;
        }
        {
            std::stringstream setname;
            setname << "xstr/" << ig;
            xsm_g.write(setname.str(), xstr, d);
        }
        {
            std::stringstream setname;
            setname << "xsnf/" << ig;
            xsm_g.write(setname.str(), xsnf, d);
        }
        {
            std::stringstream setname;
            setname << "xsch/" << ig;
            xsm_g.write(setname.str(), xsch, d);
        }
        {
            std::stringstream setname;
            setname << "xsf/" << ig;
            xsm_g.write(setname.str(), xsf, d);
        }
    }

    // scattering matrix
    VecF scat(regions_.size() * ng_ * ng_, 0.0);
    auto it = scat.begin();
    for (const auto &reg : regions_) {
        auto reg_scat = reg.xsmacsc().as_vector();
        it            = std::copy(reg_scat.begin(), reg_scat.end(), it);
    }
    VecI dims(3);
    dims[0] = regions_.size();
    dims[1] = ng_;
    dims[2] = ng_;

    xsm_g.write("xssc", scat, dims);

    return;
}
}
