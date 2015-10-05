#include "sn_boundary.hpp"

#include <iostream>

using std::cout;
using std::endl;

namespace mocc {
    SnBoundary::SnBoundary( int n_grp, int n_ang, int nx, int ny, 
            int nz ):
        n_grp_( n_grp ),
        n_ang_( n_ang ),
        nx_( nx ),
        ny_( ny ),
        nz_( nz ),
        ang_stride_( nx*ny + nx*nz + ny*nz ),
        group_stride_( ang_stride_*n_ang ),
        data_( 0.0, group_stride_* n_grp )
    {
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

    std::ostream& operator<<(std::ostream& os, const SnBoundary &b ) {
        for( size_t ig=0; ig<b.n_grp_; ig++ ) {
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
