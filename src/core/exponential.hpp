#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "core/global_config.hpp"

const int N = 10000;

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

        inline real_t exp( real_t v ) {
            return std::exp(v);
        }

        real_t max_error() {
            return max_error_;
        }

    protected:
        real_t max_error_;
    };

    /*
     * This version of \ref Exponential uses a linear lookup table to speed up
     * the evaluations of \ref exp(). If the argument to \ref exp() is beyond
     * the domain of the table, it will fall back to the result of the standard
     * library \c exp() function.
     */
    class Exponential_Linear: public Exponential {
    public:
        Exponential_Linear():
            min_( -10.0 ),
            space_( -min_/(real_t)(N-1)),
            rspace_( 1.0/space_ )
        {
            for( int i=0; i<N; i++ ) {
                d_[i] = std::exp(min_+i*space_);
            }

            for( int i=0; i<N-1; i++ ) {
                real_t x = min_+space_*(0.5+i);
                real_t err = std::abs(this->exp(x) - std::exp(x));
                max_error_ = std::max( max_error_, err );
            }
        }

        inline real_t exp( real_t v ) {
            assert((v < 0.0));
            if( v < min_ ) {
                std::cout << "Big exponential argument: " << v << std::endl;
                return std::exp(v);
            }
            int i = (v-min_)*rspace_;
            v -= space_*i+min_;
            return d_[i] + (d_[i+1] - d_[i])*v*rspace_;
        }
    protected:
        real_t min_;
        real_t space_;
        real_t rspace_;
        std::array<real_t, N> d_;
    };

    /**
     * Same as \ref Exponential_Linear, but without the check to make sure the
     * argument is within the limits of the table. This should only be used in
     * situations where one knows that the arguments to \ref exp() won't spill
     * the banks of the table. Even considering branch prediction, this manages
     * to shave a little more time off of the \ref exp() evaluations over \ref
     * Exponential_Linear.
     */
    class Exponential_UnsafeLinear : public Exponential_Linear {
        Exponential_UnsafeLinear():
            Exponential_Linear()
        {
            return;
        }

        inline real_t exp( real_t v ) {
            int i = (v-min_)*rspace_;
            v -= space_*i+min_;
            return d_[i] + (d_[i+1] - d_[i])*v*rspace_;
        }

    };
}
