#include "h5file.hpp"

#include <cassert>
#include <sstream>

#include "error.hpp"

namespace mocc {
    namespace HDF {
        H5File::H5File( std::string fname ):
            file_( fname, H5F_ACC_TRUNC )
        {
            return;
        }
    
    
        void Write( H5::CommonFG *node, std::string path, VecF data, VecI dims ) {
            hsize_t *dims_a = new hsize_t(dims.size());
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
    
    		delete dims_a;
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

        VecF::const_iterator Write( H5::CommonFG *node, std::string path,
                VecF::const_iterator first,
                VecF::const_iterator last,
                VecI dims ) 
        {
            std::vector<hsize_t> dims_a;
            int n = 1;
            for ( auto di: dims ) {
                n *= di;
                dims_a.push_back(di);
            }
            assert( (last-first) == n );

            try {
                H5::DataSpace space( dims.size(), dims_a.data() );
                H5::DataSet dataset = node->createDataSet( path, 
                        H5::PredType::NATIVE_DOUBLE, space );
                dataset.write( &(*first), H5::PredType::NATIVE_DOUBLE );
            } catch (...) {
                std::stringstream msg;
                msg << "Failed to write dataset: " << path;
                throw EXCEPT(msg.str().c_str());
            }
            return last;
        }
    }
}
