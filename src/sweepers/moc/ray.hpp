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

#include "util/global_config.hpp"
#include "geometry/geom.hpp"
#include "core_mesh.hpp"

#define MAX_NSEG 255

namespace mocc {
namespace moc {
/**
 * A \ref Ray stores vectors of segment length and the flat source region
 * index that each segment is crossing. The FSR indices are represented as
 * an offset from the first FSR in a given plane, allowing for ray data to
 * be reused for each instance of a geometrically-unique plane.
*/
class Ray {
    /**
     * This struct stores data for the "coarse ray trace," or the
     * interaction of a ray with the coarse mesh boundaries. Each entry has
     * several members, stored in a bitfield. The data on the \ref
     * RayCoarseData essentially say "move forward/backward" n segments, and
     * deposit information on the corresponding boundary. If the surface is
     * INVALID, treat the entry as a no-op.
     */
    struct RayCoarseData {
        Surface fw : 4;
        Surface bw : 4;
        unsigned int nseg_fw : 8;
        unsigned int nseg_bw : 8;

        friend std::ostream &operator<<(std::ostream &os,
                                        const RayCoarseData rcd)
        {

            os << rcd.fw << " " << rcd.nseg_fw << "\t|\t" << rcd.bw << " "
               << rcd.nseg_bw;

            return os;
        }
    };

public:
    /** \brief Construct a ray from two starting points. */
    Ray(Point2 p1, Point2 p2, std::array<int, 2> bc, int iplane,
        const CoreMesh &mesh);

    int nseg() const
    {
        return nseg_;
    }

    int ncseg() const
    {
        return cm_data_.size();
    }

    /**
     * \brief Return a reference to the coarse ray data
     */
    const std::vector<RayCoarseData> &cm_data() const
    {
        return cm_data_;
    }

    /**
     * Return the index of the first coarse mesh cell encountered by this
     * ray in the forward direction
     */
    size_t cm_cell_fw() const
    {
        return cm_cell_fw_;
    }

    /**
     * Return the index of the first coarse mesh cell encountered by this
     * ray in the backward direction.
     */
    size_t cm_cell_bw() const
    {
        return cm_cell_bw_;
    }

    /**
     * Return the index of the first coarse mesh surface encountered by this
     * ray in the forward direction
     */
    size_t cm_surf_fw() const
    {
        return cm_surf_fw_;
    }

    /**
     * Return the index of the first coarse mesh surface encountered by this
     * ray in the backward direction
     */
    size_t cm_surf_bw() const
    {
        return cm_surf_bw_;
    }

    /**
     * Return a reference to the whole vector of segment lengths
     */
    const VecF &seg_len() const
    {
        return seg_len_;
    }

    /**
     * \brief Return a reference to a single segment length.
     *
     * This method is not const, since the segment lengths must be mutable,
     * so that the \ref RayData object can correct the segment lengths once
     * all rays have been traced. Its a bit messy... when interacting with
     * the rays in any other context, the const version should be used. This
     * is relatively automatic, since the \ref RayData object only exposes
     * each \ref Ray as a const reference.
     */
    real_t &seg_len(int iseg)
    {
        return seg_len_[iseg];
    }

    /**
     * Return a const segment length
     */
    real_t seg_len(int iseg) const
    {
        return seg_len_[iseg];
    }

    /**
     * Return a reference to the whole vector of segment indices
     */
    const VecI &seg_index() const
    {
        return seg_index_;
    }

    /**
     * Only return a reference to a single segment length
     */
    size_t seg_index(size_t iseg) const
    {
        assert(iseg < nseg_);
        return seg_index_[iseg];
    }

    /**
     * Return the bc index for the start/stop of the ray
     */
    int bc(int dir) const
    {
        return bc_[dir];
    }

    /**
     * \brief Ray magnitude for > and < operators is based on number of
     * segments.
     */
    bool operator<(const Ray &other) const
    {
        return nseg_ < other.nseg_;
    }

    /**
     * \copybrief operator<()
     */
    bool operator>(const Ray &other) const
    {
        return nseg_ > other.nseg_;
    }

    friend std::ostream &operator<<(std::ostream &os, const Ray &ray);

private:
    size_t cm_surf_fw_;
    size_t cm_surf_bw_;
    size_t cm_cell_fw_;
    size_t cm_cell_bw_;

    std::vector<RayCoarseData> cm_data_;

    // Length of ray segments
    VecF seg_len_;

    // FSR index of each segment from plane offset
    VecI seg_index_;

    // Number of segments in the ray
    size_t nseg_;

    // Boundary condition index for the forward and backward directions
    std::array<int, 2> bc_;

    // Store the points that were used to initialize the ray. You know...
    // for posterity
    Point2 p1_;
    Point2 p2_;
};
}
}
