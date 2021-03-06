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
#include "pugixml.hpp"
#include "util/error.hpp"
#include "util/string_utils.hpp"
#include "pin.hpp"

using std::string;

namespace mocc {
Lattice::Lattice(const pugi::xml_node &input,
                 const std::map<int, UP_Pin_t> &pins)
{
    // Get lattice ID
    id_ = input.attribute("id").as_int(0);
    if (id_ == 0) {
        throw EXCEPT("Trouble reading lattice ID.");
    }

    // Get dimensions
    nx_ = input.attribute("nx").as_int(0);
    ny_ = input.attribute("ny").as_int(0);
    if ((nx_ == 0) | (ny_ == 0)) {
        throw EXCEPT("Trouble reading lattice dimensions.");
    }

    auto pin_ids = explode_string<int>(input.child_value());
    if (pin_ids.size() != nx_ * ny_) {
        throw EXCEPT("Incorrect number of pin IDs specified for "
                     "lattice.");
    }

    // We should be done parsing and checking things from the XML
    //
    // Arrange the pins in a 2-D array. This flips the y index from the
    // order in the input file so that the row 0, col 0 origin is in
    // the lower left.
    pins_.resize(ny_ * nx_);
    for (unsigned int iy = 0; iy < ny_; iy++) {
        unsigned int row = ny_ - iy - 1;
        for (unsigned int ix = 0; ix < nx_; ix++) {
            unsigned int col       = ix;
            int pin_id             = pin_ids[iy * nx_ + ix];
            pins_[row * nx_ + col] = pins.at(pin_id).get();
        }
    }

    // Store the pitches along each dimension
    hx_ = 0.0;
    for (unsigned int ix = 0; ix < nx_; ix++) {
        real_t dx = this->at(ix, 0).mesh().pitch_x();
        hx_ += dx;
        hx_vec_.push_back(dx);
    }

    hy_ = 0.0;
    for (unsigned int iy = 0; iy < ny_; iy++) {
        real_t dy = this->at(0, iy).mesh().pitch_y();
        hy_ += dy;
        hy_vec_.push_back(dy);
    }

    // Store the actual pin interface coordinates along each dimension
    x_vec_.push_back(0.0);
    for (unsigned int ix = 0; ix < nx_; ix++) {
        x_vec_.push_back(x_vec_[ix] + hx_vec_[ix]);
    }

    y_vec_.push_back(0.0);
    for (unsigned int iy = 0; iy < ny_; iy++) {
        y_vec_.push_back(y_vec_[iy] + hy_vec_[iy]);
    }

    // Check to make sure the pins line up nicely
    for (unsigned int iy = 0; iy < ny_; iy++) {
        for (unsigned int ix = 0; ix < nx_; ix++) {
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
    n_reg_   = 0;
    n_xsreg_ = 0;
    for (auto &pi : pins_) {
        n_reg_ += pi->mesh().n_reg();
        n_xsreg_ += pi->mesh().n_xsreg();
    }

    unsigned int prev = 0;
    first_reg_pin_.push_back(0);
    for (auto pi = pins_.begin(); pi != pins_.end() - 1; ++pi) {
        prev += (*pi)->n_reg();
        first_reg_pin_.push_back(prev);
    }
    return;
}

const PinMesh *Lattice::get_pinmesh(Point2 &p, int &first_reg,
                                    Direction dir) const
{
    assert(p.x > -REAL_FUZZ);
    assert(p.y > -REAL_FUZZ);
    assert(p.x / hx_ < 1.0 + REAL_FUZZ);
    assert(p.y / hy_ < 1.0 + REAL_FUZZ);
    // Locate the pin, and offset the point to pin-local coordinates.
    /// \todo This is potentially pretty brittle. Future PinMesh types might
    /// break the assumption here that all PinMesh origins are smack-dab in
    /// middle of the mesh. Should provide some functionality on the PinMesh
    /// itself to provide its origin to clients.
    unsigned ix = std::distance(
        x_vec_.cbegin(),
        std::lower_bound(x_vec_.cbegin(), x_vec_.cend(), p.x, fuzzy_lt));
    if (fp_equiv_abs(p.x, x_vec_[ix])) {
        ix = (dir.ox > 0.0) ? ix + 1 : ix;
    }
    ix--;
    ix = std::min(nx_ - 1, std::max(0u, ix));

    p.x = 0.5 * (x_vec_[ix + 1] + x_vec_[ix]);

    unsigned iy = std::distance(
        y_vec_.cbegin(),
        std::lower_bound(y_vec_.cbegin(), y_vec_.cend(), p.y, fuzzy_lt));
    if (fp_equiv_abs(p.y, y_vec_[iy])) {
        iy = (dir.oy > 0.0) ? iy + 1 : iy;
    }
    iy--;
    iy  = std::min(ny_ - 1, std::max(0u, iy));
    p.y = 0.5 * (y_vec_[iy + 1] + y_vec_[iy]);

    unsigned int i = iy * nx_ + ix;

    // Increment first_reg
    first_reg += first_reg_pin_[i];
    return &this->at(ix, iy).mesh();
}

std::map<int, UP_Lattice_t> ParseLattices(const pugi::xml_node &input,
                                          const std::map<int, UP_Pin_t> &pins)
{
    std::map<int, UP_Lattice_t> lattices;
    for (pugi::xml_node lat = input.child("lattice"); lat;
         lat                = lat.next_sibling("lattice")) {
        UP_Lattice_t lattice(new Lattice(lat, pins));
        if (lattices.find(lattice->id()) != lattices.end()) {
            std::stringstream msg;
            msg << "Duplicate lattice ID (" << lattice->id() << ") specified";
            throw EXCEPT(msg.str());
        }
        lattices[lattice->id()] = lattice;
    }

    return lattices;
}

bool Lattice::compatible(const Lattice &other) const
{
    if (hx_ != other.hx_) {
        std::cout << "hx " << hx_ << " " << other.hx_ << std::endl;
        return false;
    }
    if (hy_ != other.hy_) {
        std::cout << "hy" << std::endl;
        return false;
    }

    if (nx_ != other.nx_) {
        std::cout << "nx" << std::endl;
        return false;
    }
    if (ny_ != other.ny_) {
        std::cout << "ny" << std::endl;
        return false;
    }

    if (!std::equal(hx_vec_.begin(), hx_vec_.end(), other.hx_vec_.begin())) {
        std::cout << "hx_vec" << std::endl;
        return false;
    }

    if (!std::equal(hy_vec_.begin(), hy_vec_.end(), other.hy_vec_.begin())) {
        std::cout << "hy_vec" << std::endl;
        return false;
    }
    return true;
}

bool Lattice::geometrically_equivalent(const Lattice &other) const
{
    if ((nx_ != other.nx_) || (ny_ != other.ny_)) {
        return false;
    }
    if ((n_reg_ != other.n_reg_) || (n_xsreg_ != other.n_xsreg_)) {
        return false;
    }

    for(unsigned ip=0; ip<pins_.size(); ip++) {
        if(pins_[ip]->mesh().id() != other.pins_[ip]->mesh().id()) {
            return false;
        }
    }
    return true;
}

} // namespace mocc
