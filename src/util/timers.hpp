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

#include <cassert>
#include <map>

#include "global_config.hpp"

namespace mocc {
/**
 * The \ref Timer class provides functionality for measuring the amount of
 * runtime spent on various tasks. Each \ref Timer can have a number of
 * "children," which comprise sub-\ref Timer for individual tasks of
 * interest.
 *
 * Every \ref Timer maintains a total elapsed time, which may be accessed
 * via the \ref Timer::time() method. A \ref Timer can be thought of as a
 * stopwatch that is started with the \ref Timer::tic() method, and stopped
 * with the \ref Timer::toc() method. The elapsed time is a sum of all time
 * spent between calls the \ref Timer::tic() and \ref Timer::toc().
 *
 * There is a global \ref RootTimer, which is treated as the parent \ref
 * Timer for the entire executable.
 */
class Timer {
public:
    /**
     * \brief Create a new \ref Timer object.
     */
    Timer(std::string name);

    /**
     * \brief Create a new \ref Timer object and possibly start it.
     *
     * \param name the name of the \ref Timer
     * \param start whether or not to start the timer at construction-time.
     *
     * This will create a new \ref Timer and possibly start the \ref Timer
     * automatically if \p start is \c true. This is useful for
     * instrumenting constructor time for objects with lots of heavy lifting
     * to do in their initialization list. In such cases, starting the \ref
     * Timer in the body of the initializer will miss much of the time spent
     * constructing the objects in the initializer list. Placing a \ref
     * Timer object at the top of the initializer list will allow
     * measurement of this type of code.
     */
    Timer(std::string name, bool start);

    /**
     * \brief Start the \ref Timer
     *
     * This starts the \ref Timer "running," by logging the wall time at
     * which the \ref tic() function was called. The timer can then be
     * stopped with a call to \ref toc().
     */
    void tic();

    /**
     * \brief Stop the \ref Timer
     *
     * This stops the \ref Timer and returns the amount of time that elapsed
     * between the calls to \ref tic() and the \ref toc(). The running
     * sum of time for the \ref Timer is also incremented by this duration.
     */
    real_t toc();

    /**
     * \brief Return the time accumulated so far for the timer.
     *
     */
    real_t time() const;

    /**
     * \brief Return a reference the child Timer of the passed name
     */
    Timer &operator[](const std::string &name)
    {
        return children_.at(name);
    }

    /**
     * \brief Return a const reference the child Timer of the passed name
     */
    const Timer &operator[](const std::string &name) const
    {
        return children_.at(name);
    }

    /**
     * \brief Print the entire \ref Timer tree to the provided output
     * stream.
     */
    void print(std::ostream &os, int level = 0) const;

    /**
     * \brief Create a new child \ref Timer and return a reference to it.
     */
    Timer &new_timer(const std::string &name)
    {
        children_.emplace(name, Timer(name));
        return children_.at(name);
    }

    /**
     * \brief Create and return a new child \ref Timer, possibly starting it
     * automatically
     */
    Timer &new_timer(const std::string &name, bool start)
    {
        children_.emplace(name, Timer(name, start));
        return children_.at(name);
    }

    friend std::ostream &operator<<(std::ostream &os, const Timer &timer);

private:
    std::string name_;
    real_t time_;
    bool running_;
    real_t wtime_;
    std::map<std::string, Timer> children_;
};

extern Timer RootTimer;
}
