#include "h5file.hpp"

namespace mocc {
    namespace HDF {
        H5File::H5File( std::string fname ):
            file_( fname, H5F_ACC_TRUNC )
        {
            return;
        }
    
    
        void Write( H5::CommonFG *node, std::string path, VecF data, VecI dims ) {
            hsize_t *dims_a = new hsize_t(dims.size());
            for ( unsigned int i=0; i<dims.size(); i++ ) {
                dims_a[i] = dims[i];
            }
    
            H5::DataSpace space( dims.size(), dims_a );
            H5::DataSet dataset = node->createDataSet(path, 
                    H5::PredType::NATIVE_DOUBLE, space);
            dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    
    		delete dims_a;
            return;
        }
    
        void Write( H5::CommonFG *node, std::string path, int data ) {
            hsize_t dims_a[1];
    
            dims_a[0] = 1;
    
            H5::DataSpace space( 1, dims_a );
            H5::DataSet dataset = node->createDataSet(path, 
                    H5::PredType::NATIVE_INT, space);
            dataset.write(&data, H5::PredType::NATIVE_INT);
    
            return;
        }
    }
}
