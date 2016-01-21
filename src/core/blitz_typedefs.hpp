#pragma once

#include <blitz/array.h>

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

}
