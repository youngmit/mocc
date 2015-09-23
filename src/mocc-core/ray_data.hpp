#pragma once

#include <vector>
#include <memory>

#include "pugixml.hpp"

#include "core_mesh.hpp"
#include "angular_quadrature.hpp"
#include "angle.hpp"
#include "geom.hpp"

/** \file
* Here are provided a couple of classes that facilitate ray tracing for MoC.
* There is a Ray class, which contains vectors of ray segment lengths and the
* FSR index corresponding to each segment, as well as the boundary condition
* index corresponding the start and end points of the ray.
*
* Also defined is the RayData class, which is a collection of Ray objects
* organized by angle and plane.
*/


namespace mocc {
    enum VolumeCorrection {
        FLAT,
        ANGLE
    };

    /**
     * A Ray stores vectors of segment length and the flat source region index
     * that each segment is crossing. The FSR indices are represented as an
     * offset from the first FSR in a given plane, allowing for ray data to be
     * reused for each instance of a geometrically-unique plane. 
    */
    class Ray {
    public:
        /** \brief Construct a ray from two starting points. */
        Ray( Point2 p1, Point2 p2, size_t bc1, size_t bc2, int iz, 
                const CoreMesh &mesh );

        size_t nseg() const {
            return nseg_;
        }

        // Return a reference to the whole vector of segment lengths
        const VecF& seg_len() const {
            return seg_len_;
        }

        /**
         * \brief Return a reference to a single segment length.
         *
         * This method is not const, since the segment lengths must be mutable,
         * so that the RayData object can correct the segment lengths once all
         * rays have been traced. Its a bit messy... when interacting with the
         * rays in any other context, the const version should be used. This is
         * relatively automatic, since the RayData object only exposes each Ray
         * as a const reference.
         */ 
        float_t& seg_len( int iseg ) {
            return seg_len_[iseg];
        }

        /**
         * Return a const segment length
         */
        float_t seg_len( int iseg ) const {
            return seg_len_[iseg];
        }

        /**
         * Return a reference to the whole vector of segment indices
         */
        const VecI& seg_index() const {
            return seg_index_;
        }

        /**
         * Only return a reference to a single segment length
         */
        size_t seg_index( int iseg ) const {
            return seg_index_[iseg];
        }

        /**
         * Return the bc index for the start/stop of the ray
         */
        size_t bc( int dir ) const {
            return bc_[dir];
        }

        /**
         * Return the number of coarse mesh regions spanned by the ray. Corner
         * crossings technically count as two, since the current must flow
         * through adjacent cells.
         */
        size_t ncseg() const {
            return cm_cell_.size();
        }

        /** 
         * Return a reference to a vector of coarse mesh cells traversed by the
         * ray.
         *
         * \todo this needs work. Right now it doesnt conform to the surface
         * crossings
         */
        const VecI& cm_cell() const {
            assert( false );
            return cm_cell_;
        }

        /** 
         * Return the indexed coarse mesh surface crossed by the ray.
         */
        int cm_surf( int i ) const {
            return cm_surf_[i];
        }

        /** 
         * Return the number of segments on the ray 
         */
        const VecI& cm_nseg() const {
            return cm_nseg_;
        }



    private:
        // Length of ray segments
        VecF seg_len_;

        // FSR index of each segment from plane offset
        VecI seg_index_;

        // Vector containing the number of ray segments passing through each pin
        // mesh. This is needed in order to collect coarse mesh data during a
        // sweep.
        VecI cm_nseg_;

        // The coarse mesh cell index corresponding to each pin traversed by the
        // ray.
        VecI cm_cell_;

        // The coarse mesh surface index corresponding to each pin interface
        // crossed by the ray (# of pins crossed + 1).
        VecI cm_surf_;

        // Number of segments in the ray
        size_t nseg_;

        // Boundary condition index for the forward and backward directions
        size_t bc_[2];
    };

    /**
    * The RayData class is a collection of Ray objects, organized by plane, then
    * by angle.  Rays are traced only for the set of geometrically-unique planes
    * as determined by the CoreMesh object used to construct a RayData object.
    * Since the rays are only intended for use in a 2-D MoC sweeper, only the
    * first two octants are treated, with octants 3 and 4 being treated by
    * sweeping the rays backwards.
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
        typedef std::vector< std::vector <std::vector<Ray> > > RaySet_t;
    public:
        RayData( const pugi::xml_node &input, 
                 const AngularQuadrature &ang_quad,
                 const CoreMesh &mesh );

        /// Iterator to the begining of the ray data (by plane)
        RaySet_t::const_iterator begin() const {
            return rays_.cbegin();
        }

        /// Iterator to the end of the ray data (by plane)
        RaySet_t::const_iterator end() const {
            return rays_.cend();
        }

        /// Return the number of rays for the given angle index
        size_t n_rays( size_t iang ) const {
            return Nrays_[iang];
        }

        /// Return the number of rays impingent on the y-normal faces of the
        /// domain for the given angle
        size_t nx( size_t iang ) const {
            return Nx_[iang];
        }

        /// Return the number of rays impingent on the x-normal faces of the
        /// domain for the given angle
        size_t ny( size_t iang ) const {
            return Ny_[iang];
        }

        /// Return the ray spacing for the given angle
        float_t spacing( int iang ) {
            return spacing_[iang];
        }

        /// Return the maximum number of segmens spanned by any Ray in the
        /// collection. This is useful for defining the size of the scratch
        /// space for MoC.
        size_t max_segments() const {
            return max_seg_;
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
        
        // Perform a volume-correction of the ray segment lengths. This can be
        // done in two ways: using an angular integral of the ray volumes, or
        // using an angle-wice correction, which ensures that for each angle,
        // the ray segment volumes reproduce the region volumes. The first way
        // is technically more correct, however the latter is useful for
        // debugging purposes sometimes.
        void correct_volume( const CoreMesh& mesh, VolumeCorrection type );

    };

    typedef std::shared_ptr<RayData> SP_RayData_t;
}
