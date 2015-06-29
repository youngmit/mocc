#pragma once

#include <vector>

namespace mocc{

#ifdef FORCE_SINGLE
typedef float float_t;
#else
typedef double float_t;
#endif

// General purpose vector of floats
typedef std::vector<float_t> VecF;
};

#define PROG_NAME "MOCC"

