#pragma once
#include <iostream>
#include <exception>
#include <stdlib.h>
#include <string>
#include <sstream>

using std::cout;
using std::endl;
using std::string;

static void Error(const char* msg){
	cout << "ERROR: " << msg << endl;
	exit(EXIT_FAILURE);
}


static void Warn(const char* msg){
	cout << "WARNING: " << msg << endl;
}

namespace mocc {
    class Exception: public std::exception {
    public:
        Exception( const char* file, int line, const char* func, const char* msg ):
            file_( file ),
            line_( line ),
            func_( func ),
            message_( msg )
        {
            return;
        }
        const char* what() const noexcept {
            std::stringstream ret;
            ret << file_ << ":" << line_ << " in " << func_ << endl;
            ret << message_ << endl;
            return ret.str().c_str();
        }
    private:
        std::string file_;
        int line_;
        std::string func_;
        std::string message_;
    };

    static void Fail( mocc::Exception e ) {
        cout << e.what();
        exit(EXIT_FAILURE);
    }


#define EXCEPT(msg) Exception(__FILE__, __LINE__, __func__, msg);
}
