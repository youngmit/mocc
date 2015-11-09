#pragma once

#include <cassert>
#include <H5Cpp.h>
#include <memory>
#include <sstream>
#include <string>

#include "error.hpp"
#include "global_config.hpp"

namespace mocc {
    /**
     * \section hdf5_dimensions HDF5 Dataset Dimensions
     * The dimensions of the HDF5 datasets should be reversed. This is a legacy
     * carryover from my plotting software being used to interacting with data
     * from Fortran code. Make sure to flip the order of the dimensions when
     * performing output.
     */
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
         * (see \ref hdf5_dimensions )
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
         * \param dims a vector of ints containing the dimensions of the data
         * (see \ref hdf5_dimensions )
         */
        template<class InputIterator>
        InputIterator Write( H5::CommonFG *node, std::string path, 
                InputIterator first, InputIterator last, VecI dims ) {
            std::vector<hsize_t> dims_a;
            int n = 1;
            for ( auto di: dims ) {
                n *= di;
                dims_a.push_back(di);
            }

            VecF d(n);
            std::copy( first, last, d.begin() );
            assert( d.size() == n );

            try {
                H5::DataSpace space( dims.size(), dims_a.data() );
                H5::DataSet dataset = node->createDataSet( path, 
                        H5::PredType::NATIVE_DOUBLE, space );
                dataset.write( d.data(), H5::PredType::NATIVE_DOUBLE );
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
            return last;
        }

        /**
         * This is a very simple wrapper for the HDF5 file class, which
         * essentially just helps with opening a file, given a path, and can
         * return a H5::ComonFG to the root of the file.
         *
         * To write to the file (or any other reference to a H5::CommonFG) use
         * the HDF::Write functions above.
         */
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
