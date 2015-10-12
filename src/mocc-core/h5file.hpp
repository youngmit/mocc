#pragma once

#include <H5Cpp.h>
#include <memory>
#include <string>

#include "global_config.hpp"

namespace mocc {
    namespace HDF {
        /**
         * Write a vector of floats to the HDF5 file at the specified location
         * relative to the H5Group.
         *
         * \param data the vector of data. The data is read 1-dimensional, but
         * resized to the dimensions specified by \c dims.
         * \param dims a vector of ints containing the dimensions of the data
         */
        void Write( H5::CommonFG *node, std::string path, VecF data, 
                VecI dims );
        void Write( H5::CommonFG *node, std::string path, int data );

        class H5File {
        public:
            H5File( std::string fname );
            H5::CommonFG* get() {
                return &file_;
            }

        private:
            H5::H5File file_;
        };
    }
}
