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

#include "timers.hpp"

#include <iostream>
#include "util/omp_guard.h"

namespace mocc {
Timer RootTimer("MOCC");

Timer::Timer(std::string name) : name_(name), time_(0.0), running_(false)
{
    return;
}

Timer::Timer(std::string name, bool start) : Timer(name)
{
    if (start) {
        this->tic();
    }
    return;
}

void Timer::tic()
{
    assert(!running_);
    running_ = true;
    wtime_   = omp_get_wtime();
}

real_t Timer::toc()
{
    assert(running_);
    running_ = false;

    time_ += omp_get_wtime() - wtime_;
    return time_;
}

real_t Timer::time() const
{
    if (running_) {
        return time_ + omp_get_wtime() - wtime_;
    }
    else {
        return time_;
    }
}

void Timer::print(std::ostream &os, int level) const
{
    assert(!running_);
    for (int i = 0; i < level; i++) {
        os << "    ";
    }
    os << (*this) << std::endl;
    for (const auto &t : children_) {
        t.second.print(os, level + 1);
    }
    return;
}

std::ostream &operator<<(std::ostream &os, const Timer &timer)
{
    assert(!timer.running_);
    os << timer.name_ << " time: " << timer.time_ << " seconds";
    return os;
}
}
