#pragma once

#include <cassert>
#include <map>
#include <omp.h>

#include "global_config.hpp"

namespace mocc {
    /**
     * The \ref Timer class provides functionality for measuring the amount of
     * runtime spent on various tasks. Each \ref Timer can have a number of
     * "children," which comprise sub-\ref Timers for individual tasks of
     * interest.
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
        Timer(std::string name, bool start );

        void tic();

        real_t toc();

        Timer &operator[]( const std::string &name ) {
            return children_.at(name);
        }

        const Timer &operator[]( const std::string &name ) const {
            return children_.at(name);
        }

        void print( std::ostream &os, int level=0 ) const;

        Timer &new_timer( const std::string &name ) {
            children_.emplace( name, Timer(name) );
            return children_.at(name);
        }
        
        Timer &new_timer( const std::string &name, bool start ) {
            children_.emplace( name, Timer(name, start) );
            return children_.at(name);
        }

        friend std::ostream &operator<<( std::ostream &os,
                const Timer &timer );

    private:
        std::string name_;
        real_t time_;
        bool running_;
        real_t wtime_;
        std::map<std::string, Timer> children_;
    };

    extern Timer RootTimer;

}
