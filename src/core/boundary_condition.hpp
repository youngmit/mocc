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

#include <array>

#include "blitz_typedefs.hpp"
#include "angular_quadrature.hpp"
#include "constants.hpp"
#include "global_config.hpp"

namespace mocc {

    typedef std::array<int, 3> BC_Size_t;
    typedef std::array<Boundary, 6> BC_Type_t;

    typedef std::pair<int, real_t*> BVal_t;
    typedef std::pair<int, const real_t*> BVal_const_t;

    /**
     * This defines a general boundary condition container, which can store an
     * angular flux condition as a 1-dimensional collection of values for
     * each dimension-normal, angle, and energy. It is the responsibility of the
     * sweeper itself to determine what the indices of the actual conditions
     * mean.
     *
     * This class shall make one important guarantee; that all faces of the BC
     * in an angle/group are stored consecutively. This potentially allows
     * client code to eschew the concept of surface normals entirely. See \ref
     * moc::MoCSweeper for an example of this.
     */
    class BoundaryCondition {
    public:
        /**
         * \brief Construct a simple boundary condition, where there are the
         * same number of conditions per face/angle (Sn case).
         *
         * \param n_group the number of groups
         * \param angquad the \ref AngularQuadrature to use for reflections and
         * such
         * \param bc the \ref Boundary condition to enforce at each domain
         * boundary
         * \param n_bc the number of conditions to store for each direction
         * normal.
         *
         * To cut down on code duplication, just expand the scalar \p n_bc into
         * a vector and call the general case (below).
         */
        BoundaryCondition( int n_group, const AngularQuadrature &angquad,
                BC_Type_t bc, BC_Size_t n_bc );

        /**
         * \brief Construct a more complicated boundary condition, where each
         * angle can have a different number of conditions (MoC case).
         *
         * \param n_group the number of groups
         * \param angquad the \ref AngularQuadrature to use for reflections and
         * such
         * \param bc the \ref Boundary condition to enforce at each domain
         * boundary
         * \param n_bc a vector containing the number of BCs needed for each
         * angle.
         */
        BoundaryCondition( int n_group, const AngularQuadrature &angquad,
                BC_Type_t bc, std::vector< BC_Size_t > n_bc );

        /**
         * \brief Replace the default copy constructor to avoid aliasing of the
         * underlying data
         */
        BoundaryCondition( const BoundaryCondition &rhs );

        /**
         * \brief Return the total number of boundary condition points.
         */
        int size() const {
            return data_.size();
        }

        /**
         * \brief Initialize all BC points with a given value
         */
        void initialize_scalar( real_t val );

        /**
         * \brief Initialize the boundary conditions with an energy spectrum.
         */
        void initialize_spectrum( const ArrayB1 &spectrum );

        /**
         * \brief Return a const pointer to the beginning of a boundary
         * condition face
         */
        BVal_const_t get_face( int group, int angle, Normal norm ) const {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int off = bc_per_group_*group + offset_(angle, (int)norm);
            return BVal_const_t(size_[angle][(int)norm], &data_(off));
        }

        /**
         * \brief Return a pointer to the beginning of a boundary condition face
         */
        BVal_t get_face( int group, int angle, Normal norm ) {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int off = bc_per_group_*group + offset_(angle, (int)norm);
            return BVal_t(size_[angle][(int)norm], &data_(off));
        }

        /**
         * \brief Copy boundary values to external array
         */
        void copy_face( int group, int angle, Normal norm, real_t *out ) {
            int n = size_[angle][(int)norm];
            const real_t *face = this->get_face( group, angle, norm ).second;
            std::copy( face, face+n, out );
        }

        /**
         * \brief Return a const pointer to the beginning of the boundary values
         * for the given group and angle. Includes all faces
         */
        BVal_const_t get_boundary( int group, int angle ) const {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int size = size_[angle][0] + size_[angle][1] + size_[angle][2];
            int off = bc_per_group_*group + offset_(angle, 0);
            return BVal_const_t( size, &data_(off) );
        }

        /**
         * \brief Return a pointer to the beginning of the boundary values for
         * the given group and angle. Includes all faces
         */
        BVal_t get_boundary( int group, int angle ) {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int size = size_[angle][0] + size_[angle][1] + size_[angle][2];
            int off = bc_per_group_*group + offset_(angle, 0);
            return BVal_t(size, &data_(off));
        }

        /**
         * \brief Update the boundary condition for all angles for a single
         * group using a passed-in "outgoing" condition.
         *
         * \param group the energy group to update
         * \param out the "outgoing" angular flux boundary condition to use in
         * the update. It should only have one group of storage
         *
         * This would be used for a Jacobi-style iteration on the boundary
         * source.
         */
        void update( int group, const BoundaryCondition &out );

        /**
         * \brief Update the boundary condition from a single outgoing angle for
         * a single group using a passed-in "outgoing" condition.
         *
         * \param group the energy group to treat
         * \param angle the angle index of the outgoing angle. See below.
         * \param out a single-group \ref BoundaryCondition object storing the
         * outgoing boundary values to use.
         *
         * This would be used for a Gauss-Seidel-style iteration on the boundary
         * source.
         *
         * The passed angle indicates the angle of the outgoing angle. Depending
         * on the various domain boundary conditions, corresponding boundary
         * values may be updated on \c this
         */
        void update( int group, int angle, const BoundaryCondition &out );

        friend std::ostream& operator<<(std::ostream &os,
                const BoundaryCondition &bc );

    private:
        // Number of energy groups
        int n_group_;

        // Number of angles to track
        int n_angle_;

        // Boundary conditions
        std::array<Boundary, 6> bc_;

        // BC_Size_t for each angle, size_ is the same for all energy groups
        std::vector<BC_Size_t> size_;

        // Angular quadrature used to do angle index reflections
        const AngularQuadrature &ang_quad_;

        // Number of BCs per energy groups. Essentially the sum of the BCs on
        // all faces for all angles
        int bc_per_group_;

        // Large vector containing all of the boundary conditions for all
        // angles, groups and faces
        ArrayB1 data_;

        // An array of index offsets to get to an angle/face. Needs to be
        // incremented by bc_per_group_*group to yield final starting index of a
        // face of BCs
        blitz::Array<int, 2> offset_;
    };
}
