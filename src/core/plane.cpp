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
        for( auto &l: lattices_ ) {
            n_reg_   += l->n_reg();
            n_xsreg_ += l->n_xsreg();
        }
    }

    /**
     * \brief Given a Point2 in core-local coordinates, return a const pointer
     * to the corresponding PinMesh.
     * \param[in,out] p a Point2 in core-local coordinates. Will be modified
     * (see below).
     * \param[in,out] first_reg the first FSR index of the Plane. Will be
     * incremented to the first FSR index of the Lattice in which the PinMesh
     * resides.
     *
     * \note The passed Point2 \p p will be modified by the call to
     * Lattice::get_pinmesh(). See CoreMesh::get_pinmesh() for a detailed
     * description of why.
     */
    const PinMesh* Plane::get_pinmesh( Point2 &p, int &first_reg ) const {
        // Locate the lattice
        size_t ix, iy;
        for (ix=0; ix<nx_; ix++) {
            if(p.x < hx_[ix+1]) {
                break;
            }
        }
        for (iy=0; iy<ny_; iy++) {
            if(p.y < hy_[iy+1]) {
                break;
            }
        }

        size_t ilat = nx_*iy + ix;

        assert( (0<= ix) && (ix<nx_) );
        assert( (0<= iy) && (iy<ny_) );

        // Offset the point to lattice-local coordinates (distance from
        // lower-left corner of lattice)
        p.x -= hx_[ix];
        p.y -= hy_[iy];

        // Increment the first region index by the starting index of the lattice
        first_reg += first_reg_lattice_[ilat];

        // Ask lattice for reference to pin mesh, with modification of first_reg
        const PinMesh *pm = this->at(ix, iy).get_pinmesh(p, first_reg);

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
