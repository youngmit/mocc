#pragma once

#include "h5file.hpp"

namespace mocc {
    class HasOutput {
    public:
        // Perform final output to a data file. Most output is probably
        // delegated to supbordinate objects.
        virtual void output( H5File& file ) const = 0;
    };
}
