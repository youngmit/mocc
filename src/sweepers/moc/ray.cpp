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

#include "ray.hpp"

#include <cassert>

// Assuming that p1 is the "origin" return the quadrant of the angle formed by
// p1. Since we assume that p1 is below p2 in y, only octants 1 or 2 can be
// returned.
inline int get_octant(mocc::Point2 p1, mocc::Point2 p2)
{
    assert(p2.y > p1.y);
    return (p2.x > p1.x) ? 1 : 2;
}

namespace mocc {
namespace moc {
/**
 * \param p1 the starting point of the \ref Ray.
 *
 * \param p2 the ending point of the \ref Ray.
 *
 * \param bc an \c std::array<RayBC, 2> containing entries for the boundary
 * condition indices for the beginning and end of the ray.
 *
 * \param iplane the index of the geometry to trace. This corresponds to an
 * index in the \ref CoreMesh of unique geometries, which does not
 * necessarily correspond to a physical location.
 *
 * \param mesh a reference to the CoreMesh to trace.
 *
 * A Ray is defined by two \ref Point2 structs, specifying the beginning and
 * end of a ray, on the boundary of the problem. Given these two points, all
 * of the segments of the ray are determined by first finding intersections
 * with the pin cell edges (using \ref CoreMesh::trace()), then the internal
 * surface crossings for each pin (using \ref PinMesh::trace()).
 */
Ray::Ray(Point2 p1, Point2 p2, std::array<int, 2> bc, int iplane,
         const CoreMesh &mesh)
    : bc_(bc), p1_(p1), p2_(p2)
{
    std::vector<Point2> ps;

    ps.push_back(p1);
    ps.push_back(p2);

    mesh.trace(ps);
    std::vector<Point2> cps = ps;

    // Trace the fine ray. We need to keep track of the number of segments
    // in each pin crossing for the coarse ray data.
    VecI cm_nseg;
    auto p_prev = cps.front();
    for (auto pi = ps.begin() + 1; pi != ps.end(); ++pi) {
        // Use the midpoint of the pin entry and exit points to locate the
        // pin.
        auto pin_p = Midpoint(*pi, p_prev);

        int first_reg          = 0;
        // pin_p is overwritten as global coordinates of pin center.
        const PinMeshTuple pmt = mesh.get_pinmesh(pin_p, iplane, first_reg);

        int nseg = pmt.pm->trace(p_prev - pin_p, *pi - pin_p, first_reg,
                                 seg_len_, seg_index_);

        cm_nseg.push_back(nseg);

        p_prev = *pi;
    }

    // Figure out the coarse mesh data for the ray. Start with the starting
    // cells and surfaces.
    size_t ns = 0;
    Surface s[2];
    // All of the forward stuff
    int octant = get_octant(p1, p2);
    std::vector<Surface> surfs_fw;
    VecI nsegs_fw;
    {
        auto pp = cps.cbegin();
        auto p  = cps.cbegin() + 1;
        // Per convention, the first crossing is always only one surface. neat.
        cm_cell_fw_ = mesh.coarse_boundary_cell(*pp, octant);
        ns          = mesh.coarse_norm_point(*pp, octant, s);
        assert(ns == 1);
        cm_surf_fw_ = mesh.coarse_surf(cm_cell_fw_, s[0]);
        for (auto nseg : cm_nseg) {
            ns = mesh.coarse_norm_point(*p, octant, s);
            for (size_t i = 0; i < ns; i++) {
                surfs_fw.push_back(s[i]);
            }
            nsegs_fw.push_back(nseg);
            if (ns > 1) {
                nsegs_fw.push_back(0);
            }
            ++pp;
            ++p;
        }
    }
    // Backward stuff
    std::reverse(cm_nseg.begin(), cm_nseg.end());
    octant = (octant == 1) ? 3 : 4;
    std::vector<Surface> surfs_bw;
    VecI nsegs_bw;
    {
        auto pp     = cps.crbegin();
        auto p      = cps.crbegin() + 1;
        cm_cell_bw_ = mesh.coarse_boundary_cell(*pp, octant);
        ns          = mesh.coarse_norm_point(*pp, octant, s);
        assert(ns == 1);
        cm_surf_bw_ = mesh.coarse_surf(cm_cell_bw_, s[0]);
        for (auto nseg : cm_nseg) {
            ns = mesh.coarse_norm_point(*p, octant, s);
            for (size_t i = 0; i < ns; i++) {
                surfs_bw.push_back(s[i]);
            }
            nsegs_bw.push_back(nseg);
            if (ns > 1) {
                nsegs_bw.push_back(0);
            }
            ++pp;
            ++p;
        }
    }

    size_t stp = std::min(nsegs_fw.size(), nsegs_bw.size());
    for (size_t i = 0; i < stp; i++) {
        RayCoarseData rcd;
        rcd.fw      = surfs_fw[i];
        rcd.bw      = surfs_bw[i];
        rcd.nseg_fw = nsegs_fw[i];
        rcd.nseg_bw = nsegs_bw[i];

        cm_data_.push_back(rcd);
    }

    // Things get weird here. If there are different numbers of entries in
    // the forward or backward direction, ONE end of the ray must have hit a
    // corner, but not the other end. In this case, we add an extra RCD
    // instance to carry the information for the corner double-crossing
    // direction, and a no-op for the other.
    if (nsegs_fw.size() > nsegs_bw.size()) {
        RayCoarseData rcd;
        rcd.fw      = surfs_fw.back();
        rcd.bw      = Surface::INVALID;
        rcd.nseg_fw = 0;
        rcd.nseg_bw = 0;
        cm_data_.push_back(rcd);
    }
    if (nsegs_bw.size() > nsegs_fw.size()) {
        RayCoarseData rcd;
        rcd.fw      = Surface::INVALID;
        rcd.bw      = surfs_bw.back();
        rcd.nseg_fw = 0;
        rcd.nseg_bw = 0;
        cm_data_.push_back(rcd);
    }

    nseg_ = seg_len_.size();

    return;
}

std::ostream &operator<<(std::ostream &os, const Ray &ray)
{
    os << "[" << ray.p1_ << ", " << ray.p2_ << "]";

    return os;
}
}
}
