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

#include <iosfwd>
#include <memory>
#include <vector>
#include "util/pugifwd.hpp"
#include "core/angular_quadrature.hpp"
#include "core/core_mesh.hpp"
#include "core/geometry/angle.hpp"
#include "core/geometry/geom.hpp"
#include "ray.hpp"

namespace mocc {
namespace moc {
enum class VolumeCorrection { FLAT, ANGLE, NONE };
enum class Modularization { TRIG, RATIONAL };
std::ostream &operator<<(std::ostream &os, VolumeCorrection vc);
/**
 * \page coarseraypage Coarse Ray Tracing
 * Each ray crossing a mesh corner must deposit its information on one
 * exiting face of the current cell and one entering surface of the diagonal
 * neighbor. Consistency must be maintained between coincident rays of
 * different angle, otherwise surface quantities may end up with
 * nonsensical values. A good example is when current should be zero in
 * certain symmetric situations. If the corner crossings are not handled
 * properly, non-zero current could be calculated because a ray that crosses
 * one face in one direction is not being cancelled out by its sibling ray
 * in the direction reflected across that face (for instance if the
 * reflected ray passes instead through the neighboring coarse mesh
 * surface). This would impart an artificially non-zero current on both of
 * those faces.

 * \todo include discussion of coarse ray trace peculiarities and
 * conventions.
*/

/**
* The \ref RayData class is a collection of \ref Ray objects, organized by
* plane, then by angle. Rays are traced only for the set of
* geometrically-unique planes as determined by the \ref CoreMesh object used
* to construct a \ref RayData object. Since the rays are only intended for
* use in a 2-D MoC sweeper, only the first two octants are treated, with
* octants 3 and 4 being treated by sweeping the rays backwards.
*
* Boundary condition indexing is set up to be conformant with corresponding
* instances of \ref BoundaryCondition objects. The \ref BoundaryCondition
* class handles boundary values on a surface-by-surface basis, and therefore
* \ref Ray ends indexed in such a way to correspond to the appropriate faces
* on the \ref BoundaryCondition. Since the \ref BoundaryCondition stores all
* boundary values for a given angle contiguously in the X_NORM, Y_NORM,
* Z_NORM order, the ray indices should look like this:
*
* \verbatim
+- 4-- 5-- 6-- 7-- 8-- 9--10--11-+
|                                |
3                                3
|                                |
2                                2
|                                |
1                                1
|                                |
0                                0
|                                |
+- 4-- 5-- 6-- 7-- 8-- 9--10--11-+ \endverbatim
*
*/
class RayData {
    /**
     * A set of planes of \ref Ray
     */
    typedef std::vector<std::vector<Ray>> PlaneRays_t;
    typedef std::vector<PlaneRays_t> RaySet_t;

public:
    RayData(const pugi::xml_node &input, const AngularQuadrature &ang_quad,
            const CoreMesh &mesh);

    /**
     * Iterator to the beginning of the ray data (by plane)
     */
    RaySet_t::const_iterator begin() const
    {
        return rays_.cbegin();
    }

    /**
     * Iterator to the end of the ray data (by plane)
     */
    RaySet_t::const_iterator end() const
    {
        return rays_.cend();
    }

    /**
     * Return a const reference to the angular quadrature.
     *
     * The internal \ref AngularQuadrature is a modularized form of the one
     * passed in at construction time, and it is often necessary to retrieve
     * this modularized quadrature.
     */
    const AngularQuadrature &ang_quad() const
    {
        return ang_quad_;
    }

    /**
     * Return the number of rays for the given angle index
     */
    size_t n_rays(size_t iang) const
    {
        return Nrays_[iang];
    }

    /**
     * Return the number of rays impingent on the y-normal faces of the
     * domain for the given angle
     */
    size_t nx(size_t iang) const
    {
        return Nx_[iang];
    }

    /**
     * Return the number of rays impingent on the x-normal faces of the
     * domain for the given angle
     */
    size_t ny(size_t iang) const
    {
        return Ny_[iang];
    }

    /**
     * Return the ray spacing for the given angle
     */
    real_t spacing(int iang) const
    {
        return spacing_[iang];
    }

    /**
     * Return the maximum number of segments spanned by any \ref Ray in the
     * collection. This is useful for defining the size of the scratch
     * space for MoC.
     */
    int max_segments() const
    {
        return max_seg_;
    }

    /**
     * Provide stream insertion support.
     */
    friend std::ostream &operator<<(std::ostream &os, const RayData &rays);

    /**
     * \brief Return a const reference to the indexed set of plane rays.
     */
    const PlaneRays_t &operator[](size_t id) const
    {
        return rays_[id];
    }

private:
    // Methods
    std::pair<int, int> modularize_angle(Angle ang, real_t hx, real_t hy,
                                         real_t nominal_spacing) const;

    // Data
    // This starts as a copy of the angular quadrature that is passed in
    AngularQuadrature ang_quad_;

    // Vector of PlaneRays. The outer-most vector indexes the
    // geometrically-unique planes, the second index addresses the
    // individual angles, which span octants 1 and 2, and the last index
    // treats all of the rays for the given plane and angle.
    RaySet_t rays_;

    // Ray spacings for each angle. These vary from those specified due to
    // modularization
    VecF spacing_;

    // Number of rays lying on the y-normal face of the core for each angle
    VecI Nx_;
    // Number of rays lying on the x-normal face of the core for each angle
    VecI Ny_;

    // Total number of rays for a given angle
    VecI Nrays_;

    // Number of planes that we have ray data for. This is copied from
    // n_unique_planes() on the CoreMesh used to initialize the ray data.
    size_t n_planes_;

    // The type of volume correction to use
    VolumeCorrection correction_type_;

    // Maximum number of ray segments in a single ray
    int max_seg_;

    /**
     * Perform a volume-correction of the ray segment lengths. This can be
     * done in two ways: using an angular integral of the ray volumes, or
     * using an angle-wise correction, which ensures that for each angle,
     * the ray segment volumes reproduce the region volumes. The first way
     * is technically more correct, however the latter is useful for
     * debugging purposes sometimes.
     */
    void correct_volume(const CoreMesh &mesh);

    Modularization modularization_method_;
};

typedef std::shared_ptr<RayData> SP_RayData_t;
}
}
