#include "h5file.hpp"

#include <cassert>
#include <sstream>

#include "error.hpp"

using std::cout;
using std::endl;

namespace mocc {
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
