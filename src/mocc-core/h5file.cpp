#include "h5file.hpp"

namespace mocc {
    H5File::H5File( std::string fname ):
        file_(fname, H5F_ACC_TRUNC)
    {

    }

    void H5File::write( std::string path, VecF data, VecI dims ) {
        hsize_t dims_a[dims.size()];
        for ( int i=0; i<dims.size(); i++ ) {
            dims_a[i] = dims[i];
        }

        H5::DataSpace space( dims.size(), dims_a );
        H5::DataSet dataset = file_.createDataSet(path, 
                H5::PredType::NATIVE_DOUBLE, space);
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    }
}
