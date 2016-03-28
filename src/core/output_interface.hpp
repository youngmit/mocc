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
        /**
         * \brief Output relevant data to an HDF5 file node.
         */
        virtual void output( H5Node &file ) const = 0;
    };
}
