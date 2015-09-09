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

    const char* Exception::what() const noexcept {
        std::stringstream ret;
        ret << file_ << ":" << line_ << " in " << func_ << endl;
        ret << message_ << std::endl;
        return ret.str().c_str();
    }
}
