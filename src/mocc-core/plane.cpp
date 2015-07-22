#include "plane.hpp"

namespace mocc {
    Plane::Plane( const std::vector<const Lattice*> &lattices, 
            unsigned int nx, unsigned int ny ) {
        assert(lattices.size() == nx*ny);

        for (auto li=lattices.begin(); li!=lattices.end(); ++li) {
            lattices_.push_back(*li);
        }
        
        nx_ = nx;
        ny_ = ny;

        // determine the locations of the lattice interfaces
        float_t prev = 0.0;
        hx_.push_back(0.0);
        for(unsigned int ix=0; ix<nx_; ix++) {
            hx_.push_back(prev + this->at(ix, 0).hx());
            prev += this->at(ix, 0).hx();
        }
        prev = 0.0;
        hy_.push_back(0.0);
        for(unsigned int iy=0; iy<ny_; iy++) {
            hy_.push_back(prev + this->at(0, iy).hy());
            prev += this->at(0, iy).hy();
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

    const PinMesh* Plane::get_pinmesh( Point2 &p, int &first_reg ) const {
        // Locate the lattice
        unsigned int ix, iy;
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

        unsigned int ilat = nx_*iy + ix;

        assert( (0<= ix) & (ix<nx_) );
        assert( (0<= iy) & (iy<ny_) );

        // Offset the point to lattice-local coordinates (distance from
        // lower-left corner of lattice)
        p.x -= hx_[ix];
        p.y -= hy_[iy];

        // Increment the first region index by the starting index of the lattice
        first_reg += first_reg_lattice_[ilat];
        
        // Ask lattice for reference to pin mesh, with modification of first_reg
        const PinMesh* pm = this->at(ix, iy).get_pinmesh(p, first_reg);
        p.x += hx_[ix];
        p.y += hy_[iy];
        return pm;

    }
}
