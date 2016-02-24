#pragma once

#include <blitz/array.h>
#include <cassert>

#include <global_config.hpp>

namespace mocc {

#ifdef FORCE_SINGLE
typedef blitz::Array<float, 1> ArrayB1;
typedef blitz::Array<float, 2> ArrayB2;
typedef blitz::Array<float, 3> ArrayB3;
typedef blitz::Array<float, 4> ArrayB4;
#else
typedef blitz::Array<double, 1> ArrayB1;
typedef blitz::Array<double, 2> ArrayB2;
typedef blitz::Array<double, 3> ArrayB3;
typedef blitz::Array<double, 4> ArrayB4;
#endif

real_t DotProduct( const ArrayB1 &ary1, const ArrayB1 &ary2 ) {
    real_t result=0.0;
    assert( ary1.size() == ary2.size() );
    for( int i=0; i<int(ary1.size()); i++ ) {
        result = result+ary1(i)*ary2(i);
    }
    return result;
}

}
