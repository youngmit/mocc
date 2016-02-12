#pragma once

#include "blitz_typedefs.hpp"
#include "constants.hpp"
#include "error.hpp"
#include "global_config.hpp"

namespace mocc {
    
    typedef std::array<int, 3> BC_Size_t;
    typedef std::array<Boundary, 6> BC_Type_t;

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
                BC_Type_t bc, BC_Size_t n_bc ):
            BoundaryCondition( n_group,
                               angquad,
                               bc,
                               std::vector<BC_Size_t>(angquad.ndir(), n_bc)
                             )
        {
            return;
        }

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
                BC_Type_t bc, std::vector< BC_Size_t > n_bc ):
            n_group_(n_group),
            n_angle_(n_bc.size()),
            bc_(bc),
            size_(n_bc),
            ang_quad_(angquad)
        {
            assert( (angquad.ndir()   == (int)n_bc.size()) || 
                    (angquad.ndir()/2 == (int)n_bc.size()) );

            int n_angle = n_bc.size();

            bc_per_group_ = 0;
            for( auto n: n_bc ) {
                bc_per_group_ += n[0] + n[1] + n[2];
            }

            size_t total_size = bc_per_group_*n_group;
            data_.resize(total_size);
            offset_.resize(n_angle, 3);
            assert(offset_( 0, blitz::Range::all() ).isStorageContiguous());

            int offset = 0;
            int iang = 0;
            for( auto n_ang: n_bc ) {
                int iface = 0;
                for( auto n: n_ang ) {
                    offset_(iang, iface) = offset;
                    offset += n;
                    iface++;
                }
                iang++;
            }
            return;
        }

        /**
         * \brief Return the total number of boundary condition points.
         */
        int size() const {
            return data_.size();
        }

        /**
         * \brief Initialize all BC points with a given value
         */
        void initialize_scalar( real_t val ) {
            data_ = val;
        }

        /**
         * \brief Initialize the boundary conditions with an energy spectrum.
         */
        void initialize_spectrum( const ArrayB1 &spectrum ) {
            assert( (int)spectrum.size() == n_group_ );
            int it = 0;
            for( int ig=0; ig<n_group_; ig++ ) {
                real_t val = spectrum(ig);
                data_( blitz::Range(it, it+bc_per_group_-1) ) = val;
                it += bc_per_group_;
            }
            return;
        }

        /**
         * \brief Return a const pointer to the beginning of a boundary
         * condition face
         */
        const real_t* get_face( int group, int angle, Normal norm ) const {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int off = bc_per_group_*group + offset_(angle, (int)norm);
            return &data_(off);
        }

        /**
         * \brief Return a pointer to the beginning of a boundary condition face
         */
        real_t* get_face( int group, int angle, Normal norm ) {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int off = bc_per_group_*group + offset_(angle, (int)norm);
            return &data_(off);
        }

        /**
         * \brief Return a pointer to the beginning of the boundary values for
         * the given group and angle
         */
        const real_t* get_boundary( int group, int angle ) const {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int off = bc_per_group_*group + offset_(angle, 0);
            return &data_(off);
        }
        
        /**
         * \brief Return a pointer to the beginning of the boundary values for
         * the given group and angle
         */
        real_t* get_boundary( int group, int angle ) {
            assert(angle < n_angle_);
            assert(group < n_group_ );
            int off = bc_per_group_*group + offset_(angle, 0);
            return &data_(off);
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
        void update( int group, const BoundaryCondition &out ) {
            assert( out.n_group_ == 1 );

            for( int iang=0; iang<n_angle_; iang++ ) {
                this->update( group, iang, out );
            }
            return;
        }

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
        void update( int group, int angle, const BoundaryCondition &out ) {
            assert( angle < ang_quad_.ndir()/2 );
            int group_offset = bc_per_group_*group;

            for( Normal n: AllNormals ) {
                int size = size_[angle][(int)n];
                int iang_in = ang_quad_.reflect(angle, n);
                if( size == 0 ) {
                    break;
                }
                assert( size == out.size_[iang_in][(int)n]);
                assert( iang_in < n_angle_ );
                const auto &angle_in = ang_quad_[iang_in];
                int offset_in = group_offset + offset_(iang_in, (int)n);
                int offset_out = out.offset_(angle, (int)n);;

                switch( bc_[(int)(angle_in.upwind_surface(n))] ) {
                case Boundary::VACUUM:
                    data_(blitz::Range(offset_in, offset_in+size-1)) = 0.0;
                    break;

                case Boundary::REFLECT:
                    data_(blitz::Range(offset_in, offset_in+size-1)) =
                        out.data_(blitz::Range(offset_out, offset_out+size-1));
                    break;

                default:
                    throw EXCEPT("Unsupported boundary condition type");
                }
            }
            return;
        }

    private:
        // Number of energy groups
        int n_group_;

        // Number of angles to track
        int n_angle_;

        // Boundary conditions
        std::array<Boundary, 6> bc_;

        std::vector<BC_Size_t> size_;

        // Angular quadratrue used to do angle index reflections
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
