#pragma once

#include <vector>
#include <memory>

#include "pugixml.hpp"

#include "core_mesh.hpp"
#include "angular_quadrature.hpp"
#include "angle.hpp"
#include "geom.hpp"

// Here are provided a couple of classes that facilitate ray tracing for MoC.
// There is a Ray class, which contains vectors of ray segment lengths and the
// FSR index corresponding to each segment, as well as the boundary condition
// index corresponding the start and end points of the ray.
//
// The RayData class is a collection of rays, organized by plane, then by angle.
// Rays are traced only for the set of geometrically-unique planes as determined
// by the CoreMesh object used to construct a RayData object. Since the rays are
// only intended for use in a 2-D MoC sweeper, only the first two octants are
// treated, with octants 3 and 4 being treated by sweeping the rays backwards.
//
// Boundary condition indexing is somewhat arbitrary, so here's how it goes:
//
// +-17--18--19--20--21--22--23--24-+
// |                                |
// 4                                16
// |                                |
// 3                                15
// |                                |
// 2                                14
// |                                |
// 1                                13
// |                                |
// +- 5-- 6-- 7-- 8-- 9--10--11--12-+
//
// There are technically 4 angles that share a set of boundary conditions: an
// angle in quadrant 1, its reflected angle in quadrant 2, and the two angles
// pointing opposite those angles.

namespace mocc {
    class Ray {
    public:
        Ray( Point2 p1, Point2 p2, unsigned int bc1, unsigned int bc2, int iz, 
                const CoreMesh &mesh );

        unsigned int nseg() const {
            return nseg_;
        }

        // Return a reference to the whole vector of segment lengths
        const VecF& seg_len() const {
            return seg_len_;
        }
        // Only return a reference to a single segment length
        float_t& seg_len( int iseg ) {
            return seg_len_[iseg];
        }

        // Return a reference to the whole vector of segment indices
        const VecI& seg_index() const {
            return seg_index_;
        }
        // Only return a reference to a single segment length
        unsigned int seg_index( int iseg ) const {
            return seg_index_[iseg];
        }

    private:
        // Length of ray segments
        VecF seg_len_;
        // FSR index of each segment from plane offset
        VecI seg_index_;
        // Number of segments in the ray
        unsigned int nseg_;
        // Boundary condition index for the forward and backward directions
        unsigned int bc_[2];
    };


    class RayData {
        typedef std::vector< std::vector <std::vector<Ray> > > RaySet_t;
    public:
        RayData( const pugi::xml_node &input, 
                 const AngularQuadrature &ang_quad,
                 const CoreMesh &mesh );

        // Iterator to the begining of the ray data (by plane)
        RaySet_t::const_iterator begin() const {
            return rays_.cbegin();
        }

        // Iterator to the end of the ray data (by plane)
        RaySet_t::const_iterator end() const {
            return rays_.cend();
        }

        // Return the number of rays for the given angle index
        unsigned int n_rays( unsigned int iang ) const {
            return Nrays_[iang];
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
        unsigned int n_planes_;
    };

    typedef std::shared_ptr<RayData> SP_RayData_t;
}
