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

#include <cstdlib>
#include <iostream>
#include <list>
#include <sstream>

#include "error.hpp"
#include "files.hpp"

#include "util/stacktrace.hpp"

using std::cout;
using std::endl;
using std::string;

namespace mocc {

    std::list<std::string> Warnings;

    void Error(const char* msg) {
        cout << "ERROR: " << msg << endl;
        exit(EXIT_FAILURE);
    }

    void Warn(const char* msg) {
        Warnings.push_back(msg);
        LogScreen << "WARNING: " << msg << endl;
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

    Exception::Exception( const char* file, int line, const char* func,
                const std::string &msg ):
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
