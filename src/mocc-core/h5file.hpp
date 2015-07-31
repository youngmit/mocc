#pragma once

#include <string>
#include <H5Cpp.h>
#include <memory>

#include "global_config.hpp"

namespace mocc {
    // Useful typedefs
    typedef std::unique_ptr<H5::Group> UP_Group_t;


    class H5File {
    public:
        H5File(std::string fname);
        void write( std::string path, VecF data, VecI dims );
        void write( std::string path, int data );
        UP_Group_t mkdir( std::string path );
    private:
        H5::H5File file_;
    };
}
