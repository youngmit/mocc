#pragma once

#include <string>
#include <H5Cpp.h>

#include "global_config.hpp"

namespace mocc {
    class H5File {
    public:
        H5File(std::string fname);
        void write( std::string path, VecF data, VecI dims );
    private:
        H5::H5File file_;
    };
}
