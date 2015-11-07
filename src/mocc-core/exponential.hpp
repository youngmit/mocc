#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "mocc-core/global_config.hpp"

const int N = 5000;

namespace mocc {
    /**
     * Base class defining a class of simple utility classes for computing
     * exponential functions. The base version uses the stock library call to
     * exp(), while derived versions can override this with more efficient table
     * lookups.
     */
    class Exponential {
    public:
        Exponential():
            max_error_( 0.0 ) 
        { }

        inline virtual real_t exp( real_t v ) {
            return std::exp(v);
        }

        real_t max_error() {
            return max_error_;
        }

    protected:
        real_t max_error_;
    };

    class Exponential_Linear: public Exponential {
    public:
        Exponential_Linear():
            max_( -10.0 ),
            space_( max_/(real_t)N )
        {
            for( int i=0; i<N; i++ ) {
                d_[i] = std::exp(i*space_);
            }
            for( int i=0; i<N-1; i++ ) {
                real_t x = space_*(0.5+i);
                real_t err = std::abs(this->exp(x) - std::exp(x));
                max_error_ = std::max( max_error_, err );
            }
        }

        ///\todo store rspace_ instead. profile.
        inline real_t exp( real_t v ) {
            assert((v <= 0.0) && (v > max_));
            int i = v/space_;
            v -= space_*i;
            return d_[i] + (d_[i+1] - d_[i])*v/space_;
        }
    private:
        real_t max_;
        real_t space_;
        std::array<real_t, N> d_;
    };
}
