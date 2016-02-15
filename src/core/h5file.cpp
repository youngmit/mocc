#include "h5file.hpp"

#include <cassert>
#include <sstream>

#include "error.hpp"

using std::cout;
using std::endl;

unsigned int convert_access( mocc::H5Access access ) {
    switch(access) {
    case mocc::H5Access::WRITE:
        return H5F_ACC_TRUNC;
    case mocc::H5Access::APPEND:
        return H5F_ACC_RDWR;
    default:
        return H5F_ACC_RDONLY;
    }
}

namespace mocc {
    H5Node::H5Node( std::string filename, H5Access access ):
        file_( new H5::H5File( filename.c_str(), convert_access(access) ) ),
        node_( file_ ),
        access_( access )
    {
        return;
    }

    H5Node::H5Node( std::shared_ptr<H5::CommonFG> node, H5Access access ):
        file_( nullptr ),
        node_( node ),
        access_( access )
    {
        return;
    }

    H5Node H5Node::create_group( std::string path ) {
        if( access_ != H5Access::READ ) {
            H5::Group *g = new H5::Group();
            *g = node_->createGroup(path.c_str());
            std::shared_ptr<H5::CommonFG> sp( g );
            return H5Node( sp, access_ );
        } else {
            throw EXCEPT( "No write permissions" );
        }
    }

    std::vector<hsize_t> H5Node::dimensions( std::string path ) {
        try {
            H5::DataSet dataset = node_->openDataSet( path );
            H5::DataSpace dataspace = dataset.getSpace();
            size_t ndim = dataspace.getSimpleExtentNdims();
            std::vector<hsize_t> dims(ndim);
            dataspace.getSimpleExtentDims( dims.data() );

            return dims;
        }
        catch(...) {
            std::stringstream msg;
            msg << "Failed to get dataset dimensions: " << path;
            throw EXCEPT(msg.str());
        }
    }
    
    void H5Node::write( std::string path, const VecF &data ) {
        hsize_t *dims_a = new hsize_t;
        *dims_a = data.size();

        try {
            H5::DataSpace space( 1, dims_a );
            H5::DataSet dataset = node_->createDataSet( path,
                    H5::PredType::NATIVE_DOUBLE, space );
            dataset.write( data.data(), H5::PredType::NATIVE_DOUBLE );
        } catch (...) {
            std::stringstream msg;
            msg << "Failed to write dataset: " << path;
            throw EXCEPT( msg.str().c_str() );
        }

        delete dims_a;

        return;
    }

    void H5Node::write( std::string path, const VecF &data, VecI dims ) {
        hsize_t *dims_a = new hsize_t[dims.size()];
        for ( size_t i=0; i<dims.size(); i++ ) {
            dims_a[i] = dims[i];
        }

        try {
            H5::DataSpace space( dims.size(), dims_a );
            H5::DataSet dataset = node_->createDataSet( path,
                    H5::PredType::NATIVE_DOUBLE, space );
            dataset.write( data.data(), H5::PredType::NATIVE_DOUBLE );
        } catch (...) {
            std::stringstream msg;
            msg << "Failed to write dataset: " << path;
            throw EXCEPT( msg.str().c_str() );
        }

        delete[] dims_a;

        return;
    }

    void H5Node::write( std::string path, const ArrayB1 &data, VecI dims ) {
        if( !data.isStorageContiguous() ) {
            throw EXCEPT("Data is not contiguous.");
        }
        std::vector<hsize_t> dims_a;
        dims_a.reserve(dims.size());
        for( const auto &v: dims ) {
            dims_a.push_back(v);
        }

        try {
            H5::DataSpace space( dims.size(), dims_a.data() );
            H5::DataSet dataset = node_->createDataSet( path,
                    H5::PredType::NATIVE_DOUBLE, space );
            dataset.write( data.data(), H5::PredType::NATIVE_DOUBLE );
        } catch(...) {
            std::stringstream msg;
            msg << "Failed to write dataset: " << path;
            throw EXCEPT( msg.str().c_str() );
        }

        return;
    }


    void H5Node::read_1d( std::string path, ArrayB1 &data ) const {
        H5::DataSet dataset;
        H5::DataSpace dataspace;
        int h5size = -1;

        try {
            dataset = node_->openDataSet( path );
            dataspace = dataset.getSpace();
            h5size = dataspace.getSimpleExtentNpoints();
        } catch(...) {
            std::stringstream msg;
            msg << "Failed to access dataset: " << path;
            throw EXCEPT(msg.str());
        }

        if( data.size() == 0 ) {
            data.resize(h5size);
        } else {
            if( (int)data.size() != h5size ) {
                throw EXCEPT("Incompatible data sizes");
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

    void H5Node::read( std::string path, std::vector<double> &data ) const {
        H5::DataSet dataset;
        H5::DataSpace dataspace;
        int h5size = -1;
        int ndim = -1;

        try {
            dataset = node_->openDataSet( path );
            dataspace = dataset.getSpace();
            ndim = dataspace.getSimpleExtentNdims();
            h5size = dataspace.getSimpleExtentNpoints();
        } catch(...) {
            std::stringstream msg;
            msg << "Failed to access dataset: " << path;
            throw EXCEPT(msg.str());
        }

        if( ndim != 1) {
            throw EXCEPT("Vector input only supports single-dimensional data");
        }

        if( data.size() == 0 ) {
            data.resize(h5size);
        } else {
            if( (int)data.size() != h5size ) {
                throw EXCEPT("Incompatible data sizes");
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



















    namespace HDF {
        H5File::H5File( std::string fname, std::string access ) {
            if( access == "w" ) {
                file_ = H5::H5File( fname, H5F_ACC_TRUNC );
            } else if( access == "r" ) {
                file_ = H5::H5File( fname, H5F_ACC_RDONLY );
            } else {
                throw EXCEPT( "Invalid file access modality." );
            }
            return;
        }


        void Write( H5::CommonFG *node, std::string path, VecF data,
                VecI dims )
        {
            hsize_t *dims_a = new hsize_t[dims.size()];
            for ( size_t i=0; i<dims.size(); i++ ) {
                dims_a[i] = dims[i];
            }

            try {
                H5::DataSpace space( dims.size(), dims_a );
                H5::DataSet dataset = node->createDataSet(path,
                        H5::PredType::NATIVE_DOUBLE, space);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }

            delete[] dims_a;
            return;
        }

        void Write( H5::CommonFG *node, std::string path, int data ) {
            hsize_t dims_a[1];

            dims_a[0] = 1;

            try {
                H5::DataSpace space( 1, dims_a );
                H5::DataSet dataset = node->createDataSet(path,
                        H5::PredType::NATIVE_INT, space);
                dataset.write(&data, H5::PredType::NATIVE_INT);
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }

            return;
        }

        void Read( H5::CommonFG *node, std::string path, VecF &data,
                VecI &dims ) {
            try {
                H5::DataSet dataset = node->openDataSet( path );
                H5::DataSpace dataspace = dataset.getSpace();
                size_t ndim = dataspace.getSimpleExtentNdims();
                std::vector<hsize_t> my_dims(ndim);
                dataspace.getSimpleExtentDims( my_dims.data() );

                dims.clear();
                for( auto n: my_dims ) {
                    dims.push_back( n );
                }

                size_t size = 1;
                for( auto n: dims ) {
                    size *= n;
                }
                data.resize( size );

                dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE );
            } catch(...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
        }
    }
}
