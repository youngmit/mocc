#include "lattice.hpp"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "error.hpp"
#include "pin.hpp"
#include "string_utils.hpp"

using std::stringstream;
using std::string;
using std::cout;
using std::endl;

namespace mocc {
    Lattice::Lattice( const pugi::xml_node &input,
            const std::map<int, UP_Pin_t> &pins ) {
        // Get lattice ID
        id_ = input.attribute( "id" ).as_int( 0 );
        if ( id_ == 0 ) {
            Error( "Trouble reading lattice ID." );
        }

        // Get dimensions
        nx_ = input.attribute( "nx" ).as_int( 0 );
        ny_ = input.attribute( "ny" ).as_int( 0 );
        if ( (nx_ == 0) | (ny_ == 0) ) {
            Error( "Trouble reading lattice dimensions." );
        }

        // Read in the pin IDs
        {
            string pins_in = input.child_value();
            stringstream inBuf( trim(pins_in) );

            std::vector<Pin*> pin_vec;
            while ( !inBuf.eof() ) {
                int pin_id = 0;
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
                cout << pin_vec.size() << " " << nx_ << " " << ny_ << endl;
                Error( "Incorrect number of pin IDs specified for lattice." );
            }

            // We should be done parsing and checking things from the XML
            //
            // Arrange the pins in a 2-D array. This flips the y index from the
            // order in the input file so that thie row 0, col 0 origin is in
            // the lower left.
            pins_.resize( ny_*nx_ );
            for ( unsigned int iy=0; iy<ny_; iy++ ) {
                unsigned int row = ny_ - iy - 1;
                for ( unsigned int ix=0; ix<nx_; ix++ ) {
                    unsigned int col = ix;
                    pins_[row*nx_ + col] = pin_vec[iy*nx_+ix];
                }
            }

            // Store the pitches along each dimension
            hx_ = 0.0;
            for ( unsigned int ix=0; ix<nx_; ix++ ) {
                real_t dx = this->at(ix, 0).mesh().pitch_x();
                hx_ += dx;
                hx_vec_.push_back(dx);
            }

            hy_ = 0.0;
            for ( unsigned int iy=0; iy<ny_; iy++ ) {
                real_t dy = this->at(0, iy).mesh().pitch_y();
                hy_ += dy;
                hy_vec_.push_back(dy);
            }

            // Store the actual pin interface coordinates along each dimension
            x_vec_.push_back(0.0);
            for ( unsigned int ix=0; ix<nx_; ix++ ) {
                x_vec_.push_back(x_vec_[ix] + hx_vec_[ix]);
            }

            y_vec_.push_back(0.0);
            for ( unsigned int iy=0; iy<ny_; iy++ ) {
                y_vec_.push_back(y_vec_[iy] + hy_vec_[iy]);
            }

            // Check to make sure the pins line up nicely
            for ( unsigned int iy=0; iy<ny_; iy++ ) {
                for ( unsigned int ix=0; ix<nx_; ix++ ) {
                    if (this->at(ix, iy).mesh().pitch_x() != hx_vec_[ix]) {
                        Error("Inconguent pin pitches found in lattice.");
                    }
                    if (this->at(ix, iy).mesh().pitch_y() != hy_vec_[iy]) {
                        Error("Inconguent pin pitches found in lattice.");
                    }
                }
            }

            // Store the number of FSRs and XS regions
            n_reg_    = 0;
            n_xsreg_  = 0;
            for ( auto &pi: pins_ ) {
                n_reg_   += pi->mesh().n_reg();
                n_xsreg_ += pi->mesh().n_xsreg();
            }

            unsigned int prev = 0;
            first_reg_pin_.push_back(0);
            for ( auto pi=pins_.begin(); pi!=pins_.end()-1; ++pi ) {
                prev += (*pi)->n_reg();
                first_reg_pin_.push_back(prev);
            }
        }
    }

    const PinMesh* Lattice::get_pinmesh( Point2 &p, int &first_reg ) const {
        // Locate the pin, and offset the point to pin-local coordinates.
        /// \todo This is potentially pretty brittle. Future PinMesh types might
        /// breake the assumption here that all PinMesh origins are smack-dab in
        /// middle of the mesh. Should provide some functionality on the PinMesh
        /// itself to provide its origin to clients.
        unsigned int ix=0;
        unsigned int iy=0;
        for (ix=0; ix<nx_; ix++) {
            if(p.x < x_vec_[ix+1]) {
                p.x = 0.5*(x_vec_[ix+1] + x_vec_[ix]);
                break;
            }
        }
        for (iy=0; iy<ny_; iy++) {
            if(p.y < y_vec_[iy+1]) {
                p.y = 0.5*(y_vec_[iy+1] + y_vec_[iy]);
                break;
            }
        }

        unsigned int i = iy*nx_ + ix;

        // Increment first_reg
        first_reg += first_reg_pin_[i];
        return &this->at(ix, iy).mesh();
    }
}
