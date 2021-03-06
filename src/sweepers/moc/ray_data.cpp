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

#include "ray_data.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/files.hpp"
#include "util/rational_approximation.hpp"
#include "util/string_utils.hpp"
#include "util/validate_input.hpp"
#include "core/constants.hpp"
#include "core/parallel_environment.hpp"

namespace {
const std::vector<std::string> recognized_attributes = {
    "modularity", "spacing", "volume_correction", "modularization"};
}

namespace mocc {
namespace moc {
/**
 * \brief Construct a RayData object using a \<rays\> XML tag, a desired
 * AngularQuadrature, and CoreMesh.
 *
 * For now, the angular quadrature is duplicated before performing
 * modularization, which only mutates the RayData quadrature.
 *
 * \todo At some point, it is necessary to get access to the modularized
 * quadrature. There are a couple options for this:
 * - Define the modularization procedure on the AngularQuadrature itself,
 *   and perform modularization on the quadrature before we even get to this
 *   point. I don't like this, becuase it doesnt make sense conceptually to
 *   live on the quadrature itself.
 * - Provide a means to return the modularized quadrature after construction
 *   of the RayData. A little better, but perhaps cumbersome.
 * - Allow the passed AngularQuadrature reference to be non-const and mutate
 *   it directly.
 *
 * Construction performs the following steps:
 * -# Parse input from the XML
 * -# Modularize the angular quadrature and determine ray spacing for each
 *  angle
 * -# Construct Ray objects for each geometrically-unique plane and angle
 * -# Correct the ray segment lengths to preserve FSR volumes
 *
*/
RayData::RayData(const pugi::xml_node &input, const AngularQuadrature &ang_quad,
                 const CoreMesh &mesh)
    : ang_quad_(ang_quad), modularization_method_(Modularization::RATIONAL)
{
    LogScreen << "Generating ray data... " << std::endl;
    validate_input(input, recognized_attributes);

    // Make sure we have reasonable input
    if (input.empty()) {
        throw EXCEPT("No input privided for ray spacing.");
    }

    // Get the optimal ray spacing
    real_t opt_spacing = input.attribute("spacing").as_float(-1.0);
    if (opt_spacing <= 0.0) {
        throw EXCEPT("Failed to read valid ray spacing.");
    }

    // Get the volume correction option
    correction_type_ = VolumeCorrection::FLAT;
    if (!input.attribute("volume_correction").empty()) {
        std::string type = input.attribute("volume_correction").value();
        sanitize(type);
        if (type == "flat") {
            correction_type_ = VolumeCorrection::FLAT;
        } else if (type == "angle") {
            correction_type_ = VolumeCorrection::ANGLE;
        } else if (type == "none") {
            correction_type_ = VolumeCorrection::NONE;
        } else {
            throw EXCEPT("Unrecognized volume correction option in <rays>");
        }
    }

    LogScreen << "Using " << correction_type_ << " volume correction for "
                                                 "rays"
              << std::endl;

    // Get the modularity setting
    bool core_modular = true;
    if (!input.attribute("modularity").empty()) {
        std::string in_str = input.attribute("modularity").value();
        sanitize(in_str);
        if (in_str == "pin") {
            core_modular = false;
            // Make sure that all of the pins on the CoreMesh have the same
            // dimensions
            if (!mesh.is_pin_modular()) {
                throw EXCEPT(
                    "Core Mesh does not support pin modular ray "
                    "tracing.");
            }
        } else if (in_str == "core") {
            core_modular = true;
        } else {
            throw EXCEPT("Unrecognized modularity option.");
        }
    }

    if (core_modular) {
        LogFile << "Ray modularity: CORE" << std::endl;
    } else {
        LogFile << "Ray modularity: PIN" << std::endl;
    }

    // Get the modularization method
    if (!input.attribute("modularization").empty()) {
        std::string in_str = input.attribute("modularization").value();
        sanitize(in_str);
        if (in_str == "trig") {
            modularization_method_ = Modularization::TRIG;
        } else if (in_str == "rational_fraction") {
            modularization_method_ = Modularization::RATIONAL;
        } else {
            throw EXCEPT("Unrecognized modularizaion method speicified.");
        }
    }

    // Store some necessary stuff from the CoreMesh
    n_planes_ = mesh.n_unique_planes();

    // Extract whole-core dimensions
    real_t hx = mesh.hx_core();
    real_t hy = mesh.hy_core();

    // Figure out modular angles and spacings
    real_t hx_mod =
        core_modular ? mesh.hx_core() : (*mesh.begin())->mesh().pitch_x();
    real_t hy_mod =
        core_modular ? mesh.hy_core() : (*mesh.begin())->mesh().pitch_y();

    LogFile << "Original Angular quadrature " << std::endl;
    LogFile << ang_quad_ << std::endl;

    int iang = 0;
    for (auto ang_it = ang_quad_.octant(1); ang_it != ang_quad_.octant(2);
         ++ang_it) {

        Angle ang = *ang_it;

        auto N = this->modularize_angle(ang, hx_mod, hy_mod, opt_spacing);
        int Nx = N.first;
        int Ny = N.second;

        if (!core_modular) {
            Nx *= mesh.nx();
            Ny *= mesh.ny();
        }

        LogFile << "Total number of rays (Nx/Ny): " << Nx << " " << Ny
                << std::endl;

        Nx_.push_back(Nx);
        Ny_.push_back(Ny);
        Nrays_.push_back(Nx + Ny);

        real_t new_alpha = std::atan(hy * Nx / (hx * Ny));
        ang.modify_alpha(new_alpha);

        ang_quad_.modify_angle(iang, ang);
        real_t space = std::cos(ang_it->alpha) * hy / Ny;
        spacing_.push_back(space);

        iang++;
    }

    // Update weights on the angular quadrature
    ang_quad_.update_weights();

    // push more Nx, Ny, N, space onto their respective vectors, so we dont
    // have to wory about %'ing by ndir_oct
    for (iang = 0; iang < ang_quad_.ndir_oct() * 3; iang++) {
        Nx_.push_back(Nx_[iang]);
        Ny_.push_back(Ny_[iang]);
        Nrays_.push_back(Nrays_[iang]);
        spacing_.push_back(spacing_[iang]);
    }

    LogFile << "Modularized Angular quadrature " << std::endl;
    LogFile << ang_quad_ << std::endl;

    // Trace rays
    Box core_box = Box(Point2(0.0, 0.0), Point2(hx, hy));
    max_seg_     = 0;
    // loop over the planes of unique geometry
    for (unsigned iplane = 0; iplane < n_planes_; iplane++) {
        // generate rays for each angle in octants 1 and 2
        int iang = 0;
        std::vector<std::vector<Ray>> angle_rays;
        int nreg_plane = mesh.unique_plane(iplane).n_reg();
        VecI nrayfsr(nreg_plane, 0);
        for (auto ang = ang_quad_.octant(1); ang != ang_quad_.octant(3);
             ++ang) {
            int Nx  = Nx_[iang];
            int Ny  = Ny_[iang];
            std::array<int, 2> bc;
            real_t space   = spacing_[iang];
            real_t space_x = std::abs(space / std::sin(ang->alpha));
            real_t space_y = std::abs(space / std::cos(ang->alpha));

            LogFile << "Spacing: " << ang->alpha << " " << space << " "
                    << space_x << " " << space_y << std::endl;

            std::vector<Ray> rays;
            // Handle rays entering on the x-normal faces ( along the
            // y-axis)
            for (int iray = 0; iray < Ny; iray++) {
                Point2 p1;
                bc[0] = iray;
                if (ang->ox > 0.0) {
                    // We are in octant 1, enter from the left/west
                    p1.x = 0.0;
                } else {
                    // We are in octant 2, enter from the right/east
                    p1.x = hx;
                }
                p1.y      = (0.5 + iray) * space_y;
                Point2 p2 = core_box.intersect(p1, *ang);
                // The below indexing based on point position / spacing is
                // safer than it might appear at first. Since the rays are
                // laid out starting a half-spacing into the domain, the
                // i-th ray points lie between multiples of the ray spacing,
                // and are therefore sufficiently far away from multiples
                // of the spacing to permit a reliable division and cast to
                // int. Ray i will start (i+0.5)*spacing into the domain, so
                // dividing the ray position by the spacing and casting to
                // an int gives i.
                if (fp_equiv(p2.x, hx)) {
                    // BC is on the right/east boundary of the domain
                    bc[1] = p2.y / space_y;
                } else if (fp_equiv(p2.y, hy)) {
                    // BC is on the top/north boundary of the domain
                    bc[1] = p2.x / space_x + Ny;
                } else if (fp_equiv(p2.x, 0.0)) {
                    // BC is on the left/west boundary of the domain
                    bc[1] = p2.y / space_y;
                } else {
                    throw EXCEPT(
                        "Something has gone horribly wrong in the "
                        "ray trace.");
                }
                assert(bc[0] >= 0);
                assert(bc[1] >= 0);
                assert(bc[0] < Nx + Ny);
                assert(bc[1] < Nx + Ny);
                rays.emplace_back(Ray(p1, p2, bc, iplane, mesh));

                max_seg_ = std::max(rays.back().nseg(), max_seg_);
            }

            // Handle rays entering on the y-normal face
            for (int iray = 0; iray < Nx; iray++) {
                Point2 p1;
                p1.x      = (0.5 + iray) * space_x;
                p1.y      = 0.0;
                Point2 p2 = core_box.intersect(p1, *ang);
                bc[0]     = iray + Ny;
                if (fp_equiv(p2.x, hx)) {
                    // BC is on the right/east boundary of the core
                    bc[1] = p2.y / space_y;
                } else if (fp_equiv(p2.y, hy)) {
                    // BC is on the top/north boundary of the core
                    bc[1] = p2.x / space_x + Ny;
                } else if (fp_equiv(p2.x, 0.0)) {
                    // BC is on the left/west boundary of the core
                    bc[1] = p2.y / space_y;
                } else {
                    throw EXCEPT(
                        "Something has gone horribly wrong in the "
                        "ray trace.");
                }
                assert(bc[0] >= 0);
                assert(bc[1] >= 0);
                assert(bc[0] < Nx + Ny);
                assert(bc[1] < Nx + Ny);
                rays.emplace_back(Ray(p1, p2, bc, iplane, mesh));
                max_seg_ = std::max(rays.back().nseg(), max_seg_);
            }

            // Count number of ray crossings in each FSR
            for (auto &r : rays) {
                for (auto &i : r.seg_index()) {
                    nrayfsr[i]++;
                }
            }

            // Make sure that there is at least one ray in every FSR. Give a
            // warning if not.
            if (std::any_of(nrayfsr.begin(), nrayfsr.end(),
                            [](int i) { return i == 0; })) {
                Warn(
                    "No rays passed through at least one FSR. Try finer "
                    "ray spacing or larger regions.");
                for (size_t ifsr = 0; ifsr < nrayfsr.size(); ifsr++) {
                    std::cout << ifsr << " " << nrayfsr[ifsr] << std::endl;
                }
            }

            // Sort the rays by length. This might improve threading
            // performance
            // std::sort(rays.begin(), rays.end());
            // std::reverse(rays.begin(), rays.end());
            // Move the stack of rays into the vector of angular ray sets.
            angle_rays.push_back(std::move(rays));
            ++iang;
        } // Angle loop
        // Move the angular ray set to the vector of planar ray sets.
        rays_.push_back(std::move(angle_rays));
    } // Plane loop

    // Adjust ray lengths to correct FSR volume. Use an angle integral to do
    // so.
    this->correct_volume(mesh);

    LogScreen << "Done ray tracing" << std::endl;

} // RayData::RayData()

void RayData::correct_volume(const CoreMesh &mesh)
{
    switch (correction_type_) {
    // Correct each angle independently, preserving volume integral of
    // region for each angle
    case VolumeCorrection::FLAT:
        LogFile << std::endl << std::endl;
        LogFile << "Using " << correction_type_ << " volume correction for "
                << "rays." << std::endl;
        for (size_t iplane = 0; iplane < n_planes_; iplane++) {
            // flat_corr_max to store the maximum correction for all angles and
            // regions
            // flat_corr_rms to store the rms of the severity of correction
            real_t flat_corr_max = 0.0;
            int max_ireg         = 0;
            int max_iang         = 0;
            real_t flat_corr_rms = 0.0;
            const VecF &true_vol = mesh.unique_plane(iplane).areas();
            int iang             = 0;

            for (auto ang = ang_quad_.octant(1); ang != ang_quad_.octant(3);
                 ++ang) {
                VecF fsr_vol(mesh.unique_plane(iplane).n_reg(), 0.0);
                VecF flat_cf(mesh.unique_plane(iplane).n_reg(), 0.0);
                auto &rays   = rays_[iplane][iang];
                real_t space = spacing_[iang];
                for (const auto &ray : rays) {
                    for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                        size_t ireg = ray.seg_index(iseg);
                        fsr_vol[ireg] += ray.seg_len(iseg) * space;
                    }
                }

                for (size_t ireg = 0; ireg < mesh.unique_plane(iplane).n_reg();
                     ireg++) {
                    flat_cf[ireg] = true_vol[ireg] / fsr_vol[ireg];
                    // jwg
                    if (flat_corr_max < std::abs(flat_cf[ireg] - 1.0)) {
                        flat_corr_max = std::abs(flat_cf[ireg] - 1.0);
                        max_ireg      = ireg;
                        max_iang      = iang;
                    }
                    flat_corr_rms +=
                        (flat_cf[ireg] - 1.0) * (flat_cf[ireg] - 1.0);
                }

                // Correction
                for (auto &ray : rays) {
                    for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                        int ireg          = ray.seg_index(iseg);
                        ray.seg_len(iseg) = ray.seg_len(iseg) * flat_cf[ireg];
                    }
                }
                iang++;
            } // angle loop
            flat_corr_rms =
                std::sqrt(flat_corr_rms / (mesh.unique_plane(iplane).n_reg() *
                                           (ang_quad_.ndir() / 4)));

            LogFile << "For plane " << iplane
                    << ", the maximum correction occurs with "
                    << "region index " << max_ireg << " and angle index "
                    << max_iang << ", the magnitude of "
                    << "the correction being " << flat_corr_max << "."
                    << std::endl;
            LogFile << "The RMS of the correction is " << flat_corr_rms << "."
                    << std::endl;
            LogFile << std::endl << std::endl;

        } // plane loop

        break;
    // Correct all angles at the same time, preserving the angular
    // integral of the region volumes for all angles
    case VolumeCorrection::ANGLE:
        LogFile << std::endl << std::endl;
        LogFile << "Using " << correction_type_ << " volume correction for "
                << "rays." << std::endl;
        for (size_t iplane = 0; iplane < n_planes_; iplane++) {
            // flat_corr_max to store the maximum correction for all regions
            // flat_corr_rms to store the rms of the severity of correction
            real_t flat_corr_max = 0.0;
            real_t flat_corr_rms = 0.0;
            int max_ireg         = 0;

            const VecF &true_vol = mesh.unique_plane(iplane).areas();
            VecF fsr_vol(mesh.unique_plane(iplane).n_reg(), 0.0);
            int iang = 0;
            for (auto ang = ang_quad_.octant(1); ang != ang_quad_.octant(3);
                 ++ang) {
                std::vector<Ray> &rays = rays_[iplane][iang];
                real_t space           = spacing_[iang];
                real_t wgt             = ang->weight * 0.5;

                for (auto &ray : rays) {
                    for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                        int ireg = ray.seg_index(iseg);
                        fsr_vol[ireg] += ray.seg_len(iseg) * space * wgt;
                    }
                }
                ++iang;
            }
            // Convert fsr_vol into a correction factor
            for (int ireg = 0; ireg < (int)mesh.unique_plane(iplane).n_reg();
                 ireg++) {
                fsr_vol[ireg] = true_vol[ireg] / fsr_vol[ireg];
                if (flat_corr_max < std::abs(fsr_vol[ireg] - 1.0)) {
                    flat_corr_max = std::abs(fsr_vol[ireg] - 1.0);
                    max_ireg      = ireg;
                }
                flat_corr_rms += (fsr_vol[ireg] - 1.0) * (fsr_vol[ireg] - 1.0);
            }

            // Correct ray lengths to enforce proper FSR volumes
            iang = 0;
            for (auto ang = ang_quad_.octant(1); ang != ang_quad_.octant(3);
                 ++ang) {
                std::vector<Ray> &rays = rays_[iplane][iang];
                for (auto &ray : rays) {
                    for (int iseg = 0; iseg < ray.nseg(); iseg++) {
                        int ireg          = ray.seg_index(iseg);
                        ray.seg_len(iseg) = ray.seg_len(iseg) * fsr_vol[ireg];
                    }
                }
                ++iang;
            } // angle loop

            flat_corr_rms =
                sqrt(flat_corr_rms / (int)mesh.unique_plane(iplane).n_reg());

            LogFile << "For plane " << iplane
                    << ", the maximum correction occurs with "
                    << "region index " << max_ireg << ", the magnitude of "
                    << "the correction being " << flat_corr_max << "."
                    << std::endl;
            LogFile << "The RMS of the correction is " << flat_corr_rms << "."
                    << std::endl;
            LogFile << std::endl << std::endl;

        } // plane loop
        break;
    case VolumeCorrection::NONE:
        break;
    }
} // correct_volume

std::pair<int, int> RayData::modularize_angle(Angle ang, real_t hx, real_t hy,
                                              real_t nominal_spacing) const
{
    int Nx = 0;
    int Ny = 0;

    switch (modularization_method_) {
    case Modularization::RATIONAL: {
        // Get the angle cast into a slope less than 1, so that we can find a
        // rational approximation to it.
        real_t m = (ang.alpha < PI / 4.0) ? std::tan(ang.alpha)
                                          : std::tan(PI / 2.0 - ang.alpha);
        // Determine the best rational approximation of the slope.
        auto Nxy     = rational_approximation(m, 0.002, 0);
        Nx           = Nxy.first;
        Ny           = Nxy.second;
        real_t alpha = std::atan((Nx * hy) / (Ny * hx));
        int scale =
            std::ceil((hx / Nx * std::abs(std::sin(alpha))) / nominal_spacing);
        Nx *= scale;
        Ny *= scale;
        if (ang.alpha > PI / 4.0) {
            std::swap(Nx, Ny);
        }

    } break;
    case Modularization::TRIG:
    default:
        Nx = ceil(hx / nominal_spacing * std::abs(std::sin(ang.alpha)));
        Ny = ceil(hy / nominal_spacing * std::abs(std::cos(ang.alpha)));
        Nx += Nx % 2;
        Ny += Ny % 2;
        break;
    }

    return std::pair<int, int>(Nx, Ny);
}

std::ostream &operator<<(std::ostream &os, const RayData &rays)
{
    // For now, we are more interested in the rays in the macro sense. Where
    // they start and stop, more than what they do internally, so only do
    // output for one plane.

    // spit out some boilerplate python to make these easy to draw
    os << "import cairo" << std::endl;
    os << "def draw_rays( ctx, angle):" << std::endl;
    os << "    angle_rays = rays[angle]" << std::endl;

    os << "    for r in angle_rays:" << std::endl;
    os << "        p1 = r[0]" << std::endl;
    os << "        p2 = r[1]" << std::endl;
    os << "        ctx.move_to(p1[0], p1[1])" << std::endl;
    os << "        ctx.line_to(p2[0], p2[1])" << std::endl;
    os << "        ctx.close_path()" << std::endl;
    os << "    ctx.stroke()" << std::endl;
    os << "    return" << std::endl;

    const auto &plane_rays = rays.begin();
    os << "rays = [ ";
    auto angle_pos = os.tellp();
    for (auto &ang_rays : *plane_rays) {
        auto ray_pos = os.tellp();
        os << "[ " << std::endl;
        for (auto &r : ang_rays) {
            os << r;
            ray_pos = os.tellp();
            os << ",    # " << r.bc(0) << " " << r.bc(1) << std::endl;
        }
        // store the location before the comma
        auto end_pos = os.tellp();
        os.seekp(ray_pos);
        os << " ]";
        os.seekp(end_pos);
        angle_pos = os.tellp();
        os << "," << std::endl;
    }
    // go back to before the last comma and overwrite with close ]
    os.seekp(angle_pos);
    os << " ]" << std::endl;

    return os;
}

std::ostream &operator<<(std::ostream &os, VolumeCorrection vc)
{
    switch (vc) {
    case VolumeCorrection::FLAT:
        os << "FLAT";
        break;
    case VolumeCorrection::ANGLE:
        os << "ANGLE";
        break;
    case VolumeCorrection::NONE:
        os << "NONE";
        break;
    default:
        os << "UNKNOWN";
        break;
    }
    return os;
}
}
}
