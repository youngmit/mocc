#include "lattice.hpp"

#include <map>
#include <sstream>
#include <string>
#include <iostream>

#include "pin.hpp"
#include "error.hpp"
#include "string_utils.hpp"
#include "arrays.hpp"

using std::stringstream;
using std::string;

namespace mocc {
    Lattice::Lattice( const pugi::xml_node &input, 
            const std::map<int, UP_Pin_t> &pins ) {
        // Get lattice ID
        id_ = input.attribute( "id" ).as_int( -1 );
        if ( id_ == -1 ) {
            Error( "Trouble reading lattice ID." );
        }
        
        // Get dimensions
        nx_ = input.attribute( "nx" ).as_int( -1 );
        ny_ = input.attribute( "ny" ).as_int( -1 );
        if ( nx_ == -1 | ny_ == -1 ) {
            Error( "Trouble reading lattice dimensions." );
        }

        // Read in the pin IDs
        {
            string pins_in = input.child_value();
            stringstream inBuf( trim(pins_in) );

            std::vector<Pin*> pin_vec;
            while ( !inBuf.eof() ) {
                int pin_id = -1;
                inBuf >> pin_id;

                // Make sure the pin ID is valid
                if ( pins.count( pin_id ) == 0 ){
                    Error( "Unrecognized pin ID in lattice specification." );
                }
                pin_vec.push_back( pins.at( pin_id ).get() );
            }
            // Catch errors in reading
            if ( inBuf.fail() ) {
                Error( "Trouble reading pin IDs in lattice specification." );
            }

            // Make sure we have ther right number of pins
            if ( pin_vec.size() != nx_*ny_) {
                std::cout << pin_vec.size() << " " << nx_ << " " << ny_ << std::endl;
                Error( "Incorrect number of pin IDs specified for lattice." );
            }

            // We should be done parsing and checking things from the XML
            //
            // Arrange the pins in a 2-D array. This flips the y index from the
            // order in the input file so that thie row 0, col 0 origin is in
            // the lower left.
            pins_.resize( ny_*nx_ );
            for ( int iy=0; iy<ny_; iy++ ) {
                int row = ny_ - iy - 1;
                for ( int ix=0; ix<nx_; ix++ ) {
                    int col = ix;
                    pins_[row*nx_ + col] = pin_vec[iy*nx_+ix];
                }
            }

            // Store the pitches along each dimension
            hx_ = 0.0;
            for ( int ix=0; ix<nx_; ix++ ) {
                float_t dx = this->at(ix, 0)->mesh()->pitch_x();
                hx_ += dx;
                hx_vec_.push_back(dx);
            }

            hy_ = 0.0;
            for ( int iy=0; iy<ny_; iy++ ) {
                float_t dy = this->at(0, iy)->mesh()->pitch_y();
                hy_ += dy;
                hy_vec_.push_back(dy);
            }

            // Check to make sure the pins line up nicely
            for ( int iy=0; iy<ny_; iy++ ) {
                for ( int ix=0; ix<nx_; ix++ ) {
                    if (this->at(ix, iy)->mesh()->pitch_x() != hx_vec_[ix]) {
                        Error("Inconguent pin pitches found in lattice.");
                    }
                    if (this->at(ix, iy)->mesh()->pitch_y() != hy_vec_[iy]) {
                        Error("Inconguent pin pitches found in lattice.");
                    }
                }
            }
        }
    }
}
