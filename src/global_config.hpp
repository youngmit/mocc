#pragma once

#include <vector>

namespace mocc{

#ifdef FORCE_SINGLE
typedef float float_t;
#define FLOAT_EPS 1e-5
#else
typedef double float_t;
#define FLOAT_EPS 1e-12
#endif

// General purpose vector of floats
typedef std::vector<float_t> VecF;
typedef std::vector<unsigned int> VecI;
typedef std::vector<int> VecSI;
}

#define PROG_NAME "MOCC"

