#pragma once

#include <vector>
#include <memory>

#include "pugixml.hpp"

#include "core_mesh.hpp"
#include "angular_quadrature.hpp"
#include "angle.hpp"
#include "geom.hpp"

namespace mocc {
    class Ray {
    public:
        Ray( Point2 p1, Point2 p2, int iz, const CoreMesh &mesh );

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
    };


    class RayData {
    public:
        RayData( const pugi::xml_node &input, 
                 const AngularQuadrature &ang_quad,
                 const CoreMesh &mesh );
    private:
        // This starts as a copy of the angular quadrature that is passed in
        AngularQuadrature ang_quad_;

        // Vector of ray sets. The outer-most vector indexes the
        // geometrically-unique planes, the second index addresses the
        // individual angles, which span octants 1 and 2, and the last index
        // treats all of the rays for the given plane and angle.
        std::vector< std::vector< std::vector<Ray> > > rays_;

        // Ray spacings for each angle. These vary from those specified due to
        // modularization
        VecF spacing_;

        // Number of rays lying on the y-normal face of the core for each angle
        VecI Nx_;
        // Number of rays lying on the x-normal face of the core for each angle
        VecI Ny_;

        // Number of planes that we have ray data for. This is copied from
        // n_unique_planes() on the CoreMesh used to initialize the ray data.
        unsigned int n_planes_;
    };

    typedef std::shared_ptr<RayData> SP_RayData_t;
}