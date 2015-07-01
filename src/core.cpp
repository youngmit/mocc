#include "core.hpp"

#include <string>
#include <sstream>
#include <iostream>

#include "error.hpp"
#include "string_utils.hpp"

namespace mocc {
    Core::Core() {

std::cout << "default core constructor " << std::endl;
        nx_ = 0;
        ny_ = 0;
    }
    Core::Core( const pugi::xml_node &input, 
                const std::map<int, UP_Assembly_t> &assemblies):
        nx_( input.attribute("nx").as_int(-1) ),
        ny_( input.attribute("ny").as_int(-1) )
    {
std::cout << "creating core " << std::endl;
        // Make sure that we read the a proper ID
        if ( nx_ < 1 | ny_ < 1 ) {
            Error("Invalid core dimensions.");
        }

        // Read in the assembly IDs
        std::string asy_str = input.child_value();
        std::stringstream inBuf( trim(asy_str) );
        std::vector<int> asy_vec;
        for ( int i=0; i<nx_*ny_; ++i ) {
            int id;
            inBuf >> id;
            asy_vec.push_back( id );
            if (inBuf.fail()) {
                Error("Trouble reading assembly IDs in core specification.");
            }
        }
        // Store references to the assemblies in a 2D array. Make sure to flip
        // the y-index to get it into lower-left origin
        assemblies_.resize( nx_, ny_ );
        for ( int iy=0; iy<ny_; iy++ ) {
            int row = ny_ - iy - 1;
            for ( int ix=0; ix<nx_; ix++ ) {
                int col = ix;
                Assembly* asy_p = assemblies.at( asy_vec[iy*nx_+ix] ).get();
                assemblies_(row, col) = asy_p;
            }
        }
    }

    Core::~Core() {
        std::cout << "destroying core" << std::endl;
    }
}
