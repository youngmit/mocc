#include "timers.hpp"

namespace mocc {
    Timer RootTimer("MOCC");

    Timer::Timer(std::string name):
        name_(name),
        time_(0.0),
        running_(false)
    {
        return;
    }

    Timer::Timer(std::string name, bool start ):
        Timer(name)
    {
        if( start ) {
            this->tic();
        }
        return;
    }

    void Timer::tic() {
        assert(!running_);
        running_ = true;
        wtime_ = omp_get_wtime();

    }

    real_t Timer::toc() {
        assert(running_);
        running_ = false;

        time_ += omp_get_wtime() - wtime_;
        return time_;
    }

    void Timer::print( std::ostream &os, int level ) const {
        assert(!running_);
        for( int i=0; i<level; i++ ) {
            os << "    ";
        }
        os << (*this) << std::endl;
        for( const auto &t: children_ ) {
            t.second.print(os, level+1 );
        }
        return;
    }

    std::ostream &operator<<( std::ostream &os,
            const Timer &timer )
    {
        assert(!timer.running_);
        os << timer.name_ << " time: " << timer.time_ << " seconds";
        return os;
    }

    
}
