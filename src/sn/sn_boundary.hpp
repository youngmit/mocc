#pragma once

#include <cassert>
#include <iostream>

#include "constants.hpp"
#include "global_config.hpp"


namespace mocc {
    class SnBoundary {
    public:
        SnBoundary() { }
        SnBoundary( int n_grp, int n_ang , int nx, int ny, int nz);

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

        void zero_face( int grp, int ang, Normal norm ) {
            size_t start = group_stride_*grp + ang_stride_*ang;
            data_[std::slice(start, n_face_[(int)norm], 1)] = 0.0;
        }

        void initialize( float_t val ) {
            data_ = val;
        }
        
        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, 
                const SnBoundary &b );

    private:
        size_t n_grp_;
        size_t n_ang_;
        size_t nx_;
        size_t ny_;
        size_t nz_;
        int ang_stride_;
        int group_stride_;
        int face_offset_[3];
        int n_face_[3];
        ArrayF data_;
    };
}
