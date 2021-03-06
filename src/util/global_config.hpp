/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <vector>

namespace mocc {

#ifdef FORCE_SINGLE
typedef float real_t;
typedef int32_t int_t;
#else
typedef double real_t;
typedef int64_t int_t;
#endif

// General purpose vector of floats, ints, etc
typedef std::vector<real_t> VecF;
typedef std::vector<int> VecI;
typedef std::vector<int> VecSI;
}

#define PROG_NAME "MOCC"
