#include <cassert>
#include <iostream>
#include <vector>

namespace mocc {
    template <class T> class Array2D {
    public:
        Array2D() {
            d1_ = 0;
            d2_ = 0;
        }

        Array2D( size_t d1, size_t d2 ):
            d1_(d1), d2_(d2)
        {
            assert(d1_ * d2_ > 0);
            data_ = std::vector<T>(d1_*d2_);
        }

        //
        Array2D(const Array2D &x) {
            d1_   = x.d1();
            d2_   = x.d2();
            data_ = x.data();
            for (size_t i=0; i<d1_*d2_; i++) {
                data_[i] = x(i);
            }
        }

        ~Array2D() {
        }

        size_t d1() const {
            return d1_;
        }

        size_t d2() const {
            return d2_;
        }

        std::vector<T>& data() const {
            return data_;
        }   

        Array2D& resize( size_t new_d1 , size_t new_d2 ) {
            assert( new_d1 > 0 & new_d2 > 0 );
            size_t n = new_d1*new_d2;
            if ( n != d1_*d2_ ) {
                data_.resize( n );
                d1_ = new_d1;
                d2_ = new_d2;
            }
            return *this;
        }

        // flattened subscript operator
        T& operator()( size_t i ) {
            return data_[i];
        }

        // flattened subscript operator (const)
        const T& operator()( size_t i ) const {
            return data_[i];
        }

        // Full-dimension subscript operation
        T& operator()( size_t i, size_t j ) {
            return data_[d2_*j + i];
        }
        
        // Full-dimension subscript operation (const)
        const T& operator()( size_t i, size_t j ) const {
            return data_[d2_*j + i];
        }


        // Assignment (deep copy)
        Array2D<T>& operator=( const Array2D<T> &rhs ) {
            // Check for self-assignment
            if ( &rhs == this ){
                return *this;
            }

            d1_ = rhs.d1();
            d2_ = rhs.d2();

            data_ = rhs.data();
            return *this;
        }

        // Scalar assignment
        Array2D<T>& operator=( T &val ) {
            for( auto &v: data_ ) {
                v = val;
            }
            return *this;
        }

        typename std::vector<T>::iterator row_begin( size_t row ) {
            return data_.begin() + row*d2_;
        }

        typename std::vector<T>::iterator row_end( size_t row ) {
            return data_.begin() + (row+1)*d2_;
        }

        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, 
                const Array2D &array) {
            for( size_t j=0; j<array.d2(); j++ ) {
                for( size_t i=0; i<array.d1(); i++ ) {
                    os << array(i, j) << " ";   
                }
                os << std::endl;
            }
            return os;
        }

    private:
        size_t d1_;
        size_t d2_;
        std::vector<T> data_;
    };
}
