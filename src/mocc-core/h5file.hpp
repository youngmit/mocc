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
         * \param node the HDF5 file 
         * \param path the path to the dataset, relative to \c node. If
         * preceeded by a '/', the path is absolute and will resolve to a
         * location relative to the root of the HDF5 file.
         * \param data the vector of data. The data is read 1-dimensional, but
         * resized to the dimensions specified by \c dims.
         * \param dims a vector of ints containing the dimensions of the data
         */
        void Write( H5::CommonFG *node, std::string path, VecF data, 
                VecI dims );

        /**
         * Write a single integer to the HDF5 file at the specified location
         * relative to the H5Group.
         *
         * \param node the location in the file to write to
         * \param path the path to the dataset, relative to \c node. If
         * preceeded by a '/', the path is absolute and will resolve to a
         * location relative to the root of the HDF5 file.
         * \param data the int to be written.
         */
        void Write( H5::CommonFG *node, std::string path, int data );

        /**
         * Write data to and HDF5 location using iterators.
         *
         * \param node the location in the file to write to
         * \param path the path to the dataset, relative to \c node. If
         * preceeded by a '/', the path is absolute and will resolve to a
         * location relative to the root of the HDF5 file.
         * \param first an iterator pointing to the beginning of the data to be
         * written.
         * \param last an iterator past the end of the data to be written.
         */
        VecF::const_iterator Write( H5::CommonFG *node, std::string path, 
                VecF::const_iterator first,
                VecF::const_iterator last,
                VecI dims );

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
