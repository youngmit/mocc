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

#include "parallel_environment.hpp"
#include <omp.h>
#include "pugixml.hpp"
#include "util/error.hpp"

namespace mocc {
ParallelEnvironment::ParallelEnvironment(const pugi::xml_node &input)
    : num_threads_(1)
{
    if (!input.empty()) {
        num_threads_ =
            input.attribute("num_threads").as_int(0);

        if (num_threads_ < 1) {
            throw EXCEPT(
                "Less than one thread specified in <parallel> "
                "tag");
        }

        if (num_threads_ > omp_get_num_procs()) {
            Warn(
                "More threads specified than physical "
                "threads on this machine in <parallel> tag");
        }

        // Okay. Looking good. Tell OpenMP whats up
        omp_set_num_threads(num_threads_);
    }
    return;
}

ParallelEnvironment ParEnv;
}
