#include "sn_boundary.hpp"

#include <iostream>

using std::cout;
using std::endl;

namespace mocc {
    SnBoundary::SnBoundary( int n_group, const AngularQuadrature &ang_quad,
            const Mesh &mesh ):
        n_group_( n_group ),
        ang_quad_( ang_quad ),
        n_ang_( ang_quad.ndir() ),
        nx_( mesh.nx() ),
        ny_( mesh.ny() ),
        nz_( mesh.nz() ),
        ang_stride_( nx_*ny_ + nx_*nz_ + ny_*nz_ ),
        group_stride_( ang_stride_*n_ang_ ),
        bc_( mesh.boundary() ),
        data_( 0.0, group_stride_* n_group )
    {
        assert( n_group_ > 0 );
        n_face_[(int)Normal::X_NORM] = ny_*nz_;
        n_face_[(int)Normal::Y_NORM] = nx_*nz_;
        n_face_[(int)Normal::Z_NORM] = nx_*ny_;
        int offset = 0;
        for( const auto norm: AllNormals ) {
            face_offset_[(int)norm] = offset;
            offset += n_face_[(int)norm];
        }

        return;
    }

    void SnBoundary::update( size_t group, const SnBoundary &out ) {
        assert( group < n_group_ );
        assert( out.n_group() == 1 );

        // In the below, iang is the incoming angle index to be set, iang_refl
        // is the reflected angle from which a reflective BC should get its flux
        int iang = 0;
        for( auto &ang: ang_quad_ ) {
            Surface surf;
            Normal norm;
            int iang_refl;

            // X-normal surface
            norm = Normal::X_NORM;
            surf = (ang.ox > 0) ? Surface::WEST : Surface::EAST;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_[(int)surf] == Boundary::REFLECT ) {
                auto new_bc = out.get_face( 0, iang_refl, norm );
                this->set_face( group, iang, norm, new_bc );
            } else {
                this->zero_face(group, iang, norm);
            }
            // Y-normal surface
            norm = Normal::Y_NORM;
            surf = (ang.oy > 0) ? Surface::SOUTH : Surface::NORTH;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_[(int)surf] == Boundary::REFLECT ) {
                auto new_bc = out.get_face( 0, iang_refl, norm );
                this->set_face(group, iang, norm, new_bc);
            } else {
                this->zero_face(group, iang, norm);
            }
            // Z-normal surface
            norm = Normal::Z_NORM;
            surf = (ang.oz > 0) ? Surface::BOTTOM : Surface::TOP;
            iang_refl = ang_quad_.reflect( iang, norm );
            if( bc_[(int)surf] == Boundary::REFLECT ) {
                auto new_bc = out.get_face( 0, iang_refl, norm );
                this->set_face(group, iang, norm, new_bc );
            } else {
                this->zero_face(group, iang, norm);
            }

            iang++;
        }
        return;
    }

    void SnBoundary::update( size_t group, size_t ang, const SnBoundary &out ) {
        Angle angle = ang_quad_[ang];
        Surface surf;
        Normal norm;
        int iang_refl = 0;

        /// \todo replace with a loop over AllNormal
        // X-normal surface
        norm = Normal::X_NORM;
        surf = (angle.ox > 0) ? Surface::WEST : Surface::EAST;
        iang_refl = ang_quad_.reflect( ang, norm );
        if( bc_[(int)surf] == Boundary::REFLECT ) {
            auto new_bc = out.get_face( 0, ang, norm );
            this->set_face( group, iang_refl, norm, new_bc );
        } else {
            this->zero_face(group, iang_refl, norm);
        }
        // Y-normal surface
        norm = Normal::Y_NORM;
        surf = (angle.oy > 0) ? Surface::SOUTH : Surface::NORTH;
        iang_refl = ang_quad_.reflect( ang, norm );
        if( bc_[(int)surf] == Boundary::REFLECT ) {
            auto new_bc = out.get_face( 0, ang, norm );
            this->set_face(group, iang_refl, norm, new_bc);
        } else {
            this->zero_face(group, iang_refl, norm);
        }
        // Z-normal surface
        norm = Normal::Z_NORM;
        surf = (angle.oz > 0) ? Surface::BOTTOM : Surface::TOP;
        iang_refl = ang_quad_.reflect( ang, norm );
        if( bc_[(int)surf] == Boundary::REFLECT ) {
            auto new_bc = out.get_face( 0, ang, norm );
            this->set_face(group, iang_refl, norm, new_bc );
        } else {
            this->zero_face(group, iang_refl, norm);
        }
    }

    std::ostream& operator<<(std::ostream& os, const SnBoundary &b ) {
        for( size_t ig=0; ig<b.n_group_; ig++ ) {
            cout << "Group " << ig << endl;
            for( size_t iang=0; iang<b.n_ang_; iang++ ) {
                cout << "Angle " << iang << endl;
                cout << "X-Normal face:" << endl;
                ArrayF f = b.get_face( ig, iang, Normal::X_NORM );
                for( auto v: f ) {
                    cout << v << " ";
                }
                cout << endl;
                cout << "Y-Normal face:" << endl;
                f = b.get_face( ig, iang, Normal::Y_NORM );
                for( auto v: f ) {
                    cout << v << " ";
                }
                cout << endl;
                cout << "Z-Normal face:" << endl;
                f = b.get_face( ig, iang, Normal::Z_NORM );
                for( auto v: f ) {
                    cout << v << " ";
                }
                cout << endl;
                cout << endl;

            }
        }
        return os;
    }
}
