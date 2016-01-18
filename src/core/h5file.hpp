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

        H5Node create_group( std::string path );

        H5::CommonFG *get() {
            return node_.get();
        }

        H5Node operator[]( std::string path ) {
            std::shared_ptr<H5::CommonFG> 
                g(new H5::Group(node_->openGroup(path)));
            return H5Node( g, access_ );
        }

        void write( std::string path, VecF data, VecI dims );

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

        template<class BlitzArray>
        void read( std::string path, BlitzArray &data ) {
            try {
                H5::DataSet dataset = node_->openDataSet( path );
                H5::DataSpace dataspace = dataset.getSpace();
                size_t ndim = dataspace.getSimpleExtentNdims();
                std::vector<hsize_t> dims(ndim);
                dataspace.getSimpleExtentDims( dims.data() );
                
                if( (int)ndim != data.dimensions() ) {
                    throw EXCEPT("Array and dataset dimensionality disagree.");
                }

                auto shape = data.shape();
                auto size = data.size();

                if( size == 0 ) {
                    // Allocate space in the destination array
                    for( unsigned i=0; i<ndim; i++ ) {
                        shape[i] = dims[i];
                    }
                    data.resize(shape);
                } else {
                    // Make sure the dimensions match
                    for( unsigned i=0; i<ndim; i++ ) {
                        if( shape[i] != (int)dims[i] ) {
                            throw EXCEPT("Incorrect array shape.");
                        }
                    }
                }

                dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE );
            } catch(...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
            return;
        }

        /**
         * Write data to an HDF5 location using iterators.
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
        InputIterator write( H5::CommonFG *node, std::string path,
                InputIterator first, InputIterator last, VecI dims ) {
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
         * Write data to an HDF5 location using iterators.
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
            assert( (int)d.size() == n );

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


        void Read( H5::CommonFG *node, std::string path, VecF &data,
                VecI &dims );

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
            /**
             * \brief Open a new HDF5 file.
             *
             * \param fname the path to the file to open
             * \param access the access modality to use.
             *
             * Current options for \p access  are "r" and "w". If using "r," the
             * file is opened as-is and no writes can be performed. When using
             * "w" any existing file is erased and replaced with the data
             * written.
             */
            H5File( std::string fname, std::string access );
            H5::CommonFG* get() {
                return &file_;
            }

        private:
            H5::H5File file_;
        };
    }
}
