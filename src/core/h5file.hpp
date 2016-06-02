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

#include <cassert>
#include <H5Cpp.h>
#include <memory>
#include <sstream>
#include <string>

#include "error.hpp"
#include "blitz_typedefs.hpp"
#include "global_config.hpp"

namespace mocc {
    /**
     * Enumeration of the supported HDF5 access patterns supported by \ref
     * H5Node.
     */
    enum class H5Access {
        /**
         * Leave file as-is with read-only permissions.
         */
        READ,
        /**
         * Delete any existing file, opening a new one with write permissions.
         */
        WRITE,
        /**
         * Open a file as-is, with read-write permissions.
         */
        APPEND,
    };

    /**
     * Wrapper class for the HDF5 library. An instance of this class is similar
     * to the HDF5 CommonFG class, with extra stuff to make it easier to work
     * with.
     */
    class H5Node {
    public:
        H5Node( std::string filename, H5Access access );

        /**
         * \brief Create a new group in the HDF5 file.
         *
         * \param path the path to the group to be created
         * \return a new \ref H5Node pointing to the new group
         *
         * This method produces a new HDF5 group in the file, located at the
         * path specified relative to the \ref H5Node upon which this method is
         * called. If the path contains a leading slash ('/'), the path will be
         * relative to the root of the file.
         */
        H5Node create_group( std::string path );

        /**
         * \brief Return a pointer to the underlying H5::CommonFG
         *
         * This is provided so that client code can use the full-blown HDF5
         * library to interact with the file.
         */
        H5::CommonFG *get() {
            return node_.get();
        }

        /**
         * \brief Return an \ref H5Node pointing to the path specified.
         */
        H5Node operator[]( std::string path ) {
            std::shared_ptr<H5::CommonFG>
                g(new H5::Group(node_->openGroup(path)));
            return H5Node( g, access_ );
        }

        /**
         * \brief Return the dimensions of the dataset specified by the path
         */
        std::vector<hsize_t> dimensions( std::string path );

        /**
         * \brief Write an std::vector<double> to the file, assuming
         * one dimension.
         */
        void write( std::string path, const VecF &data );

        /**
         * \brief Write an std::vector<double> to the file, using specified
         * dimensions.
         */
        void write( std::string path, const VecF &data, VecI dims );

        void write( std::string path, const std::string &str );

        /**
         * \brief Write a 1-D Blitz array to an HDF5 file, possibly reshaping
         */
        void write( std::string path, const ArrayB1 &data, VecI dims );

        /**
         * \brief Write a Blitz++ array to the HDF5 node.
         *
         * This can be used to write a Blitz array of any dimension to an HDF5
         * file. This passes the pointer to the beginning of the data and
         * allows the HDF5 library to just copy directly, so it is required that
         * the data be stored contiguously, and will throw an exception if this
         * isnt the case. We also assume that the data type stored in the Blitz
         * array is double precision.
         */
        template<class BlitzArray>
        void write( std::string path, const BlitzArray &data ) {
            // First, make sure that the data is contiguous in memory. We are
            // going to do a direct copy, using the data() pointer to the
            // beginning of the array, so this is crucial.
            if( !data.isStorageContiguous() ) {
                throw EXCEPT("Blitz data is not contiguous.");
            }

            int ndim = data.dimensions();
            std::vector<hsize_t> dims;
            size_t size = 1;
            for( int i=0; i<ndim; i++ ) {
                dims.push_back( data.extent(i) );
                size *= data.extent(i);
            }

            if( data.size() != size ) {
                throw EXCEPT("Something fishy is going on with array "
                        "dimensions.");
            }

            try {
                H5::DataSpace space( dims.size(), dims.data() );
                H5::DataSet dataset = node_->createDataSet(path,
                        H5::PredType::NATIVE_DOUBLE, space);
                dataset.write( data.data(), H5::PredType::NATIVE_DOUBLE );
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }

            return;
        }

        /**
         * \brief Version of \ref H5Node::write() for writing scalar integers.
         */
        void write( std::string path, int data ) {
            hsize_t dims_a[1];

            dims_a[0] = 1;

            try {
                H5::DataSpace space( 1, dims_a );
                H5::DataSet dataset = node_->createDataSet(path,
                        H5::PredType::NATIVE_INT, space);
                dataset.write(&data, H5::PredType::NATIVE_INT);
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
            return;
        }

        /**
         * \brief Write a scalar unsigned long integer.
         */
        void write( std::string path, unsigned long data ) {
            hsize_t dims_a[1];

            dims_a[0] = 1;

            try {
                H5::DataSpace space( 1, dims_a );
                H5::DataSet dataset = node_->createDataSet(path,
                        H5::PredType::NATIVE_ULONG, space);
                dataset.write(&data, H5::PredType::NATIVE_ULONG);
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
            return;
        }

        /**
         * \brief Read data from an \ref H5Node into an STL vector
         */
        void read( std::string path, std::vector<double> &data ) const;

        /**
         * \brief Read data from an \ref H5Node into a 1-D Blitz array
         *
         * Dimensionality of the data in the HDF5 file is allowed to be more
         * than 1-D, in which case the data is copied linearly into the 1-D
         * Blitz array. This is not an overload of the regular \c read() method
         * to avoid the potential for error.
         */
        void read_1d( std::string path, ArrayB1 &data ) const;

        /**
         * \brief Read data from an \ref H5Node into a Blitz++ array.
         *
         *
         * This will attempt to read the dataset specified by the path relative
         * the location of the \ref H5Node into the passed Blitz array.
         * \pre One of the following must be true:
         *  - The array is one-dimensional, or
         *  - The array has the same dimensions as the HDF5 database.
         * \pre One of the following musht be true:
         *  - The array is empty (size == 0)
         *  - The array has space allocated for the same number of elements as
         *  the HDF5 dataset, and \c isStorageContiguous() must return \c true.
         */
        template<class BlitzArray>
        void read( std::string path, BlitzArray &data ) const {
            if( (data.size() > 0) && !data.isStorageContiguous() ) {
                throw EXCEPT("Blitz data is not contiguous");
            }

            H5::DataSet dataset;
            H5::DataSpace dataspace;
            int ndim = -1;
            int h5size = -1;
            std::vector<hsize_t> dims;
            try {
                dataset = node_->openDataSet( path );
                dataspace = dataset.getSpace();
                ndim = dataspace.getSimpleExtentNdims();
                h5size = dataspace.getSimpleExtentNpoints();
                dims.resize(ndim);
                dataspace.getSimpleExtentDims( dims.data() );
            } catch(...) {
                std::stringstream msg;
                msg << "Failed to access dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }

            if( (int)ndim != data.dimensions() ) {
            }

            auto shape = data.shape();
            auto size = data.size();

            if( size == 0 ) {
                // Allocate space in the destination array
                if( data.dimensions() == 1) {
                    data.resize(h5size);
                } else {
                    if(data.dimensions() != ndim ) {
                        throw EXCEPT("Array and dataset dimensionality "
                                "disagree.");
                    }
                    for( int i=0; i<ndim; i++ ) {
                        shape[i] = dims[i];
                    }
                    data.resize(shape);
                }

            } else {
                if( data.dimensions() == 1 ) {
                    if( (int)data.size() != h5size ) {
                        std::cerr << data.size() << " " << h5size << std::endl;
                        throw EXCEPT("Incorrect array size.");
                    }
                } else {
                    // Make sure the dimensions match
                    for( int i=0; i<ndim; i++ ) {
                        if( shape[i] != (int)dims[i] ) {
                            throw EXCEPT("Incorrect array shape.");
                        }
                    }
                }
            }

            try {
                dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE );
            } catch(...) {
                std::stringstream msg;
                msg << "Failed to read dataset: " << path;
                throw EXCEPT(msg.str());
            }

            return;
        }

        /**
         * Write data to an HDF5 location using iterators.
         *
         * \param path the path to the dataset, relative to \c node. If
         * preceeded by a '/', the path is absolute and will resolve to a
         * location relative to the root of the HDF5 file.
         * \param first an iterator pointing to the beginning of the data to be
         * written.
         * \param last an iterator past the end of the data to be written.
         * \param dims a vector of ints containing the dimensions of the data
         */
        template<class InputIterator>
        InputIterator write( std::string path, InputIterator first,
                InputIterator last, VecI dims )
        {
            std::vector<hsize_t> dims_a;
            int n = 1;
            for( auto di: dims ) {
                n *= di;
                dims_a.push_back(di);
            }

            VecF d(n);
            std::copy( first, last, d.begin() );
            assert( (int)d.size() == n );

            try {
                H5::DataSpace space( dims.size(), dims_a.data() );
                H5::DataSet dataset = node_->createDataSet( path,
                        H5::PredType::NATIVE_DOUBLE, space );
                dataset.write( d.data(), H5::PredType::NATIVE_DOUBLE );
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
            return last;
        }

    private:
        H5Node( std::shared_ptr<H5::CommonFG> node_, H5Access access );

        // Pointer to the file object. Null if not the root node of the file.
        std::shared_ptr< H5::H5File > file_;
        std::shared_ptr< H5::CommonFG > node_;
        H5Access access_;
    };
}
