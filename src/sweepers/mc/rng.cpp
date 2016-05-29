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

#include "rng.hpp"

#include <omp.h>

namespace mocc {
namespace mc {
/**
 * This is the executable-global random number generator, used by any of the
 * monte carlo components that need it. If ever there were a reason to have
 * global data, this is it.
 */
RNGSwarm RNG_SWARM(0, omp_get_max_threads());

} // namespace mc
} // namespace mocc
