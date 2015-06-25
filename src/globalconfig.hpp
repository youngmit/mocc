#pragma once

namespace mocc{

#ifdef FORCE_SINGLE
typedef float float_t;
#else
typedef double float_t;
#endif
};

#define PROG_NAME "MOCC"