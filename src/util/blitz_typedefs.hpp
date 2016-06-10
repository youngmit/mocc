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

#include <blitz/array.h>
#include <cassert>
#include "global_config.hpp"

namespace mocc {

typedef blitz::Array<real_t, 1> ArrayB1;
typedef blitz::Array<real_t, 2> ArrayB2;
typedef blitz::Array<real_t, 3> ArrayB3;
typedef blitz::Array<real_t, 4> ArrayB4;

inline real_t DotProduct(const ArrayB1 &ary1, const ArrayB1 &ary2)
{
    real_t result = 0.0;
    assert(ary1.size() == ary2.size());
    for (int i = 0; i < int(ary1.size()); i++) {
        result = result + ary1(i) * ary2(i);
    }
    return result;
}
}
