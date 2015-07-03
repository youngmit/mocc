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
        Ray( Point2 p, Angle ang );

        int nseg() const {
            return nseg_;
        }
    private:
        // Length of ray segments
        VecF seg_len_;
        // FSR index of each segment from plane offset
        VecI seg_index_;
        // Number of segments in the ray
        int nseg_;
    };


    class RayData {
    public:
        RayData( const pugi::xml_node &input, 
                 const AngularQuadrature &ang_quad,
                 const CoreMesh &mesh );
    private:
        // This starts as a copy of the angular quadrature that is passed in
        // to the constructor, but gets mutated to 
        AngularQuadrature ang_quad_;

        // Vector of ray objects
        std::vector<Ray> rays_;

        // Ray spacings for each angle. These vary from those specified due to
        // modularization
        VecF spacing_;
    };

    typedef std::shared_ptr<RayData> SP_RayData_t;
}
