#pragma once

#include <vector>
#include <memory>

#include "pugixml.hpp"

#include "angle.hpp"
#include "angular_quadrature.hpp"
#include "core_mesh.hpp"
#include "geom.hpp"
#include "ray.hpp"

namespace mocc {
    enum VolumeCorrection {
        FLAT,
        ANGLE
    };

    /**
    * The \ref RayData class is a collection of \ref Ray objects, organized by
    * plane, then by angle. Rays are traced only for the set of
    * geometrically-unique planes as determined by the \ref CoreMesh object used
    * to construct a \ref RayData object.  Since the rays are only intended for
    * use in a 2-D MoC sweeper, only the first two octants are treated, with
    * octants 3 and 4 being treated by sweeping the rays backwards.
    *
    * Boundary condition indexing is somewhat arbitrary, so here's how it goes:
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
    * There are technically 4 angles that share a set of boundary conditions: an
    * angle in quadrant 1, its reflected angle in quadrant 2, and the two angles
    * pointing opposite those angles.
    */
    class RayData {
        /**
         * A set of planes of \ref Ray
         */
        typedef std::vector< std::vector <std::vector<Ray> > > RaySet_t;
        typedef std::vector< std::vector<Ray> > PlaneRays_t;
    public:
        RayData( const pugi::xml_node &input,
                 const AngularQuadrature &ang_quad,
                 const CoreMesh &mesh );

        /**
         * Iterator to the begining of the ray data (by plane)
         */
        RaySet_t::const_iterator begin() const {
            return rays_.cbegin();
        }

        /**
         * Iterator to the end of the ray data (by plane)
         */
        RaySet_t::const_iterator end() const {
            return rays_.cend();
        }

        /**
         * Return a const reference to the angular quadratre.
         *
         * The internal \ref AngularQuadrature is a modularized form of the one
         * passed in at construction time, and it is often necessary to retrieve
         * this modularized quadrature.
         */
        const AngularQuadrature &ang_quad() const {
            return ang_quad_;
        }

        /**
         * Return the number of rays for the given angle index
         */
        size_t n_rays( size_t iang ) const {
            return Nrays_[iang];
        }

        /**
         * Return the number of rays impingent on the y-normal faces of the
         * domain for the given angle
         */
        size_t nx( size_t iang ) const {
            return Nx_[iang];
        }

        /**
         * Return the number of rays impingent on the x-normal faces of the
         * domain for the given angle
         */
        size_t ny( size_t iang ) const {
            return Ny_[iang];
        }

        /**
         * Return the ray spacing for the given angle
         */
        real_t spacing( int iang ) const {
            return spacing_[iang];
        }

        /**
         * Return the maximum number of segments spanned by any \ref Ray in the
         * collection. This is useful for defining the size of the scratch
         * space for MoC.
         */
        size_t max_segments() const {
            return max_seg_;
        }

        /**
         * Provide stream insertion support.
         */
        friend std::ostream& operator<<( std::ostream &os,
                const RayData &rays );

        /**
         * \brief Return a const reference to the indexed set of plane rays.
         */
        const PlaneRays_t& operator[]( size_t id ) const {
            return rays_[id];
        }

    private:
        // This starts as a copy of the angular quadrature that is passed in
        AngularQuadrature ang_quad_;

        // Vector of ray sets. The outer-most vector indexes the
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

        // Maximum number of ray segments in a single ray
        size_t max_seg_;

        /**
         * Perform a volume-correction of the ray segment lengths. This can be
         * done in two ways: using an angular integral of the ray volumes, or
         * using an angle-wice correction, which ensures that for each angle,
         * the ray segment volumes reproduce the region volumes. The first way
         * is technically more correct, however the latter is useful for
         * debugging purposes sometimes.
         */
        void correct_volume( const CoreMesh& mesh, VolumeCorrection type );
    };

    typedef std::shared_ptr<RayData> SP_RayData_t;
}
