/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "lattice.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>

#include "error.hpp"
#include "pin.hpp"
#include "string_utils.hpp"

using std::string;
using std::cout;
using std::endl;

namespace mocc {
    Lattice::Lattice( const pugi::xml_node &input,
            const std::map<int, UP_Pin_t> &pins ) {
        // Get lattice ID
        id_ = input.attribute( "id" ).as_int( 0 );
        if ( id_ == 0 ) {
            throw EXCEPT( "Trouble reading lattice ID." );
        }

        // Get dimensions
        nx_ = input.attribute( "nx" ).as_int( 0 );
        ny_ = input.attribute( "ny" ).as_int( 0 );
        if ( (nx_ == 0) | (ny_ == 0) ) {
            throw EXCEPT( "Trouble reading lattice dimensions." );
        }

        auto pin_ids = explode_string<int>(input.child_value());
        if( pin_ids.size() != nx_*ny_ ) {
            throw EXCEPT( "Incorrect number of pin IDs specified for "
                    "lattice." );
        }

        // We should be done parsing and checking things from the XML
        //
        // Arrange the pins in a 2-D array. This flips the y index from the
        // order in the input file so that the row 0, col 0 origin is in
        // the lower left.
        pins_.resize( ny_*nx_ );
        for ( unsigned int iy=0; iy<ny_; iy++ ) {
            unsigned int row = ny_ - iy - 1;
            for ( unsigned int ix=0; ix<nx_; ix++ ) {
                unsigned int col = ix;
                int pin_id = pin_ids[iy*nx_+ix];
                pins_[row*nx_ + col] = pins.at( pin_id ).get();
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
                    throw EXCEPT("Incongruent pin pitches found in "
                            "lattice.");
                }
                if (this->at(ix, iy).mesh().pitch_y() != hy_vec_[iy]) {
                    throw EXCEPT("Incongruent pin pitches found in "
                            "lattice.");
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
        return;
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

    std::map<int, UP_Lattice_t> ParseLattices( const pugi::xml_node &input,
            const std::map<int, UP_Pin_t> &pins ) {
        std::map<int, UP_Lattice_t> lattices;
        for ( pugi::xml_node lat = input.child( "lattice" ); lat;
                lat = lat.next_sibling( "lattice" )) {
            UP_Lattice_t lattice( new Lattice(lat, pins) );
            lattices[lattice->id()] = lattice;
        }

        return lattices;
    }

    bool Lattice::compatible( const Lattice &other ) const {
        if( hx_ != other.hx_ ) {
            cout << "hx " << hx_ << " " << other.hx_ << endl;
            return false;
        }
        if( hy_ != other.hy_ ) {
            cout << "hy" << endl;
            return false;
        }

        if( nx_ != other.nx_ ) {
            cout << "nx" << endl;
            return false;
        }
        if( ny_ != other.ny_ ) {
            cout << "ny" << endl;
            return false;
        }

        if( !std::equal( hx_vec_.begin(), hx_vec_.end(),
                    other.hx_vec_.begin() ) )
        {
            cout << "hx_vec" << endl;
            return false;
        }

        if( !std::equal( hy_vec_.begin(), hy_vec_.end(),
                    other.hy_vec_.begin() ) )
        {
            cout << "hy_vec" << endl;
            return false;
        }
        return true;
    }
}
