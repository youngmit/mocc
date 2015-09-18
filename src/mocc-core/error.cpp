#include <sstream>

#include "error.hpp"

using std::cout;
using std::endl;
using std::string;

namespace mocc {
    void Error(const char* msg) {
	    cout << "ERROR: " << msg << endl;
	    exit(EXIT_FAILURE);
    }


    void Warn(const char* msg) {
    	cout << "WARNING: " << msg << endl;
    }
    
    void Fail( Exception e ) {
        std::cout << e.what();
        exit(EXIT_FAILURE);
    }

    Exception::Exception( const char* file, int line, const char* func, 
                const char* msg ):
        file_( file ),
        line_( line ),
        func_( func ),
        message_( msg )
    {
        std::stringstream ret;
        ret << file_ << ":" << line_ << " in " << func_ << endl;
        ret << message_ << std::endl;
        print_message_ = ret.str();

        return;
    }


    const char* Exception::what() const noexcept {
        return print_message_.c_str();
    }
}
