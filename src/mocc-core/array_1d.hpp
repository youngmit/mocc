namespace mocc {
    template <class T> class Array1D {
    public:
        Array1D() {
            d1_ = 0;
        }

        Array1D( size_t d1 ):
            d1_(d1)
        {
            assert(d1_ > 0);
            data_ = std::vector<T>(d1_);
        }

        //
        Array1D(const Array1D &x) {
            d1_   = x.d1();
            data_ = x.data();
            for (size_t i=0; i<d1_; i++) {
                data_[i] = x(i);
            }
        }

        ~Array1D() {
        }

        size_t d1() const {
            return d1_;
        }

        std::vector<T>& data() const {
            return data_;
        }

        Array2D& resize( size_t new_d1  ) {
            assert( new_d1 > 0 );
            size_t n = new_d1;
            if ( n != d1_ ) {
                data_.resize( n );
                d1_ = new_d1;
            }
            return *this;
        }

        // flattened subscript operator
        T& operator()( size_t i ) const {
            return data_[i];
        }

        // Assignment (deep copy)
        Array1D<T>& operator=( const Array1D<T> &rhs ) {
            // Check for self-assignment
            if ( &rhs == this ){
                return *this;
            }

            d1_ = rhs.d1();

            data_ = rhs.data();
            return *this;
        }

        // Scalar assignment
        Array1D<T>& operator=( T &val ) {
            for( auto &v: data_ ) {
                v = val;
            }
            return *this;
        }

    private:
        size_t d1_;
        std::vector<T> data_;
    };
}
