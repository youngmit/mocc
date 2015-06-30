#pragma once 

#include <cassert>
#include <iostream>

namespace mocc {
    template <class T> class Array2D {
    public:
        Array2D() {
            d1_ = 0;
            d2_ = 0;
        }

        Array2D( int d1, int d2 ):
            d1_(d1), d2_(d2)
        {
            data_ = new T[d1_*d2_];
        }

        Array2D(const Array2D &x) {
            d1_ = x.d1();
            d2_ = x.d2();
            data_ = new T[d1_*d2_];
            for (int i=0; i<d1_*d2_; i++) {
                data_[i] = x(i);
            }
        }

        ~Array2D() {
            if ( d1_ > 0 & d2_ > 0 ) {
                delete[] data_;
            }
            d1_ = 0;
            d2_ = 0;
        }

        int d1() const {
            return d1_;
        }

        int d2() const {
            return d2_;
        }

        Array2D<T>& resize( int new_d1 , int new_d2 ) {
            int n = new_d1*new_d2;
            if ( n != d1_*d2_ ) {
                if ( d1_*d2_ > 0 ) {
                    delete[] data_;
                }
                data_ = new T[n];
                d1_ = new_d1;
                d2_ = new_d2;
            }
            return *this;
        }

        // flattened subscript operator
        T& operator()( int i ) const {
            return data_[i];
        }

        // Full-dimension subscript operation
        T& operator()( int i, int j ) const {
            return data_[d2_*i+j];
        }

        // Assignment (deep copy)
        Array2D<T>& operator=( const Array2D<T> &rhs ) {
            // Check for self-assignment
            if ( &rhs == this ){
                return *this;
            }

            // See if we need to reallocate
            int n = rhs.d1()*rhs.d2();
            if ( d1_*d2_ != n ) {
                if( d1_*d2_ > 0 ) {
                    delete[] data_;
                }
                data_ = new T[n];
            }

            d1_ = rhs.d1();
            d2_ = rhs.d2();

            for( int i=0; i<n; i++ ) {
                data_[i] = rhs(i);
            }

            return *this;
        }

        

    private:
        int d1_;
        int d2_;
        T* data_;
    };
}
