#pragma once
#include <iostream>
#include <exception>
#include <string>

namespace mocc {
    extern void Error(const char* msg);

    extern void Warn(const char* msg);



    class Exception: public std::exception {
    public:
        Exception( const char* file, int line, const char* func,
                const char* msg );

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
