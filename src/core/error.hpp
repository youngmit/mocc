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
#include <exception>
#include <string>

namespace mocc {
    extern void Error(const char* msg);

    extern void Warn(const char* msg);



    class Exception: public std::exception {
    public:
        Exception( const char* file, int line, const char* func,
                const char* msg );
        Exception( const char* file, int line, const char* func,
                const std::string &msg );

        const char* what() const noexcept;

    private:
        std::string file_;
        int line_;
        std::string func_;
        std::string message_;
        std::string print_message_;
    };


    extern void Fail( Exception e );


#define EXCEPT(msg) Exception(__FILE__, __LINE__, __func__, msg);
}
