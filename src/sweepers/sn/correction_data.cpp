#include "correction_data.hpp"

#include <iomanip>

using std::setfill;
using std::setw;

namespace mocc {
    void CorrectionData::output( H5::CommonFG *file ) const {
        VecI dims;
        dims.push_back(nz_);
        dims.push_back(ny_);
        dims.push_back(nx_);

        file->createGroup("/alpha_x");
        file->createGroup("/alpha_y");
        file->createGroup("/beta");

        for( size_t g=0; g<ngroup_; g++ ) {
            for( size_t a=0; a<nang_; a++ ) {
                VecF alpha_x( nreg_, 0.0 );
                VecF alpha_y( nreg_, 0.0 );
                VecF beta( nreg_, 0.0 );
                for( size_t i=0; i<nreg_; i++ ) {
                    alpha_x[i] = this->alpha( i, a, g,
                            Normal::X_NORM );
                    alpha_y[i] = this->alpha( i, a, g,
                            Normal::Y_NORM );
                    beta[i] = this->beta( i, a, g );
                }

                {
                    std::stringstream setname;
                    setname << "/beta/" << g << "_" << a;
                    HDF::Write( file, setname.str(), beta, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_x/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                    HDF::Write( file, setname.str(), alpha_x, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_y/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                    HDF::Write( file, setname.str(), alpha_y, dims );
                }
            }
        }
        return;
    }
    
}
