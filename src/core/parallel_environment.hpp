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

#include "util/pugifwd.hpp"

namespace mocc {
/**
 * \brief Small class for storing information about the parallel environment
 * in which the code is running.
 *
 * There should be a global instance of this class, initialized somewhere
 * like \ref InputProc or similar. For now, since we only do threading this
 * class is not much more useful than a global integer storing the number of
 * threads available. As we get more complicated or add MPI, this will
 * become more useful.
 *
 * Really, this should be a singleton class, but I have better things to do
 * than make sure I handle all the nasty corner cases that that would
 * entail.
 */
class ParallelEnvironment {
private:
    int num_threads_;

public:
    ParallelEnvironment() : num_threads_(1)
    {
        return;
    }

    ParallelEnvironment( const pugi::xml_node &input);

    int num_threads() const
    {
        return num_threads_;
    }

    void set_num_threads(int num_threads)
    {
        num_threads_ = num_threads;
    }
};

// Declare the global instance of ParallelEnvironment
extern ParallelEnvironment ParEnv;
}
