#pragma once

#include "core_mesh.hpp"
#include "geom.hpp"
#include "global_config.hpp"

namespace mocc {
    /**
     * A Ray stores vectors of segment length and the flat source region index
     * that each segment is crossing. The FSR indices are represented as an
     * offset from the first FSR in a given plane, allowing for ray data to be
     * reused for each instance of a geometrically-unique plane. 
    */
    class Ray {
        struct RayCoarseData {
            unsigned int i; 
        };
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
        size_t cm_start_fw_;
        size_t cm_start_bw_;
        
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
}
