#pragma once

#include "h5file.hpp"

namespace mocc {
    /**
     * This simply specifies an interface for outputing "stuff" to and HDF5
     * file. Any class extending it must implement the output() method, which
     * adds its data to the H5File instance passed to it.
     */
    class HasOutput {
    public:
        // Perform final output to a data file. Most output is probably
        // delegated to supbordinate objects.
        virtual void output( H5::CommonFG *file ) const = 0;
    };
}
