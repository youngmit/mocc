#pragma once

#include <cassert>
#include <iostream>

#include "core/angular_quadrature.hpp"
#include "core/constants.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"


namespace mocc {
    class SnBoundary {
    public:
        SnBoundary( int n_group, const AngularQuadrature &ang_quad,
                const Mesh &mesh );

        ArrayF get_face( int grp, int ang, Normal norm ) const {
            int in = (int)norm;
            size_t start = group_stride_*grp + ang_stride_*ang +
                face_offset_[in];
            return data_[std::slice( start, n_face_[in], 1)];
        }

        void set_face( int grp, int ang, Normal norm, const ArrayF &in ) {
            size_t start = group_stride_*grp + ang_stride_*ang +
                face_offset_[(int)norm];
            size_t size = n_face_[(int)norm];
            assert(in.size() == size);
            data_[std::slice(start, size, 1)] = in;
        }

        /**
         * \brief Apply a zero boundary condition to an entire face.
         */
        void zero_face( int grp, int ang, Normal norm ) {
            size_t start = group_stride_*grp + ang_stride_*ang +
                face_offset_[(int)norm];
            size_t size = n_face_[(int)norm];
            data_[std::slice(start, size, 1)] = 0.0;
        }

        /**
         * \brief Initialize the boundary condition to a single value.
         */
        void initialize( real_t val ) {
            data_ = val;
        }

        /**
         * \brief Return the number of energy groups for which the \ref
         * SnBoundary is defined.
         */
        size_t n_group() const {
            return n_group_;
        }

        /**
         * Update an incoming boundary condition for all angles of a given group
         *
         * \param[in] group the group to update
         * \param[in] out the outgoing angular flux boundary condition to use
         * for the update. Should only be allocated to a single group.
         *
         * This will perform a loop over all angles in the \ref SnBoundary and
         * use the passed-in "outgoing" BC to update the state of the "incoming"
         * BC. This method should be called on the incoming \ref SnBoundary
         * object and should be passed the outgoing condition. It is assumed
         * that the outgoing condition is only allocated to a single group.
         *
         * This is useful for BC updates when performing a Jacobi-like iteration
         * in the angle space.
         */
        void update( size_t group, const SnBoundary &out );

        /**
         * Update an incoming boundary condition for a single angle of a given
         * group
         *
         * \param[in] group the group to update
         *
         * \param[in] ang the angle to update from. For non-vacuum conditions,
         * BCs are transfered TO the reflection of \p ang FROM \p ang.
         *
         * \param[in] out the outgoing angular flux boundary condition to use
         * for the update. Should only be allocated to a single group.
         *
         * This will update the state of the "incoming" BC. This method should
         * be called on the incoming \ref SnBoundary object and should be passed
         * the outgoing condition. It is assumed that the outgoing condition is
         * only allocated to a single group.
         *
         * This is useful for BC updates when performing a Gauss-Seidel-lide
         * iteration in the angle space.
         */
        void update( size_t group, size_t ang, const SnBoundary &out );


        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os,
                const SnBoundary &b );

    private:
        size_t n_group_;
        AngularQuadrature ang_quad_;
        size_t n_ang_;
        size_t nx_;
        size_t ny_;
        size_t nz_;
        int ang_stride_;
        int group_stride_;
        int face_offset_[3];
        int n_face_[3];
        std::vector<Boundary> bc_;
        ArrayF data_;
    };
}
