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

#include "plane.hpp"

#include <iostream>

#include "error.hpp"
#include "fp_utils.hpp"

namespace mocc {
    Plane::Plane( const std::vector<const Lattice*> &lattices,
            size_t nx, size_t ny ):
        lattices_(lattices)
    {
        assert(lattices.size() == nx*ny);

        nx_ = nx;
        ny_ = ny;

        // determine the locations of the lattice interfaces
        VecF dx;
        real_t prev = 0.0;
        hx_.push_back(0.0);
        for( size_t ix=0; ix<nx_; ix++) {
            hx_.push_back(prev + this->at(ix, 0).hx());
            dx.push_back(this->at(ix, 0).hx());
            prev += this->at(ix, 0).hx();
        }
        VecF dy;
        prev = 0.0;
        hy_.push_back(0.0);
        for( size_t iy=0; iy<ny_; iy++ ) {
            hy_.push_back(prev + this->at(0, iy).hy());
            dy.push_back(this->at(0, iy).hy());
            prev += this->at(0, iy).hy();
        }

        // Ensure that the lattices in the plane conform
        for( size_t iy=0; iy<ny_; iy++ ) {
            for( size_t ix=0; ix<nx_; ix++ ) {
                if( !fp_equiv_ulp( this->at(ix, iy).hx(), dx[ix] )) {
                    throw EXCEPT("Lattices do not have compatible dimensions.");
                }
                if( !fp_equiv_ulp( this->at(ix, iy).hy(), dy[iy] )) {
                    throw EXCEPT("Lattices do not have compatible dimensions.");
                }
            }
        }

        // Store the first region index for each lattice within the plane
        first_reg_lattice_.push_back(0);
        auto prev_reg = 0;
        for (auto li=lattices_.begin(); li!=lattices_.end()-1; ++li) {
            prev_reg += (*li)->n_reg();
            first_reg_lattice_.push_back(prev_reg);
        }

        // Accumulate the number of FSRs and XS mesh regions
        n_reg_   = 0;
        n_xsreg_ = 0;
        n_fuel_  = 0;
        for( auto &l: lattices_ ) {
            n_reg_   += l->n_reg();
            n_xsreg_ += l->n_xsreg();
            for( const auto pin: *l ) {
                if( pin->is_fuel() ) {
                    n_fuel_++;
                }
            }
        }

        return;
    }

    /**
     * \note The passed Point2 \p p will be modified by the call to
     * Lattice::get_pinmesh(). See CoreMesh::get_pinmesh() for a detailed
     * description of why.
     */
    const PinMesh* Plane::get_pinmesh( Point2 &p, int &first_reg,
            Direction dir) const {
        assert(p.x > -REAL_FUZZ);
        assert(p.y > -REAL_FUZZ);
        assert(p.x/hx_.back() < 1.0+REAL_FUZZ);
        assert(p.y/hy_.back() < 1.0+REAL_FUZZ);
        // Locate the lattice
        unsigned ix = std::distance(
            hx_.begin(),
            std::lower_bound(hx_.begin(), hx_.end(), p.x, fuzzy_lt));
        if(fp_equiv_abs(p.x, hx_[ix])) {
            ix = (dir.ox > 0.0) ? ix+1 : ix;
        }
        ix--;
        unsigned iy = std::distance(
            hy_.begin(),
            std::lower_bound(hy_.begin(), hy_.end(), p.y, fuzzy_lt));
        if(fp_equiv_abs(p.y, hy_[iy])) {
            iy = (dir.oy > 0.0) ? iy+1 : iy;
        }
        iy--;

        // Force indices to be inside the plane
        ix = std::min(nx_-1, std::max(0u, ix));
        iy = std::min(ny_-1, std::max(0u, iy));

        size_t ilat = nx_*iy + ix;

        // Offset the point to lattice-local coordinates (distance from
        // lower-left corner of lattice)
        p.x -= hx_[ix];
        p.y -= hy_[iy];

        // Increment the first region index by the starting index of the lattice
        first_reg += first_reg_lattice_[ilat];

        // Ask lattice for reference to pin mesh, with modification of first_reg
        const PinMesh *pm = this->at(ix, iy).get_pinmesh(p, first_reg, dir);

        // Restore the point coordinates to core-local
        p.x += hx_[ix];
        p.y += hy_[iy];
        return pm;
    }

    Position Plane::pin_position( size_t ipin ) const {
        int ilat = 0;
        for( auto &lattice: lattices_ ) {
            if( ipin < lattice->n_pin() ) {
                break;
            }
            ipin -= lattice->n_pin();
            ilat++;
        }


        // ilat should be the index of the lattice in which the pin resides.
        // ipin should be the lattice-local index of the pin
        Position pos(0, 0, 0);

        int lat_x = ilat % nx_;
        int lat_y = ilat / nx_;

        for( int ix=0; ix<lat_x; ix++ ) {
            pos.x += this->at(ix, 0).nx();
        }
        pos.x += ipin % this->at(lat_x, lat_y).nx();

        for( int iy=0; iy<lat_y; iy++ ) {
            pos.y += this->at(0, iy).ny();
        }
        pos.y += ipin / this->at(lat_x, lat_y).nx();

        return pos;
    }
}
