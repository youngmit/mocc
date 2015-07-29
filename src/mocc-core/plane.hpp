#pragma once

#include <vector>
#include <iostream>

#include "lattice.hpp"
#include "geom.hpp"
#include "global_config.hpp"

namespace mocc {
    class Plane {
    public:
        Plane(const std::vector<const Lattice*> &lattices, unsigned int nx, 
                unsigned int ny);
        
        const Lattice& at(unsigned int ix, unsigned int iy) const {
            return *(lattices_[ix + nx_*iy]);
        }

        const PinMesh* get_pinmesh( Point2 &p, int &first_reg) const;        

        unsigned int n_reg() const {
            return n_reg_;
        }

        unsigned int n_xsreg() const {
            return n_xsreg_;
        }

        const VecF vol() const {
            VecF vol;
            for( auto &lat: lattices_ ) {
                for( auto &pin: *lat ) {
                    vol.insert(vol.end(), pin->vol().begin(), pin->vol().end());
                }
            }

            return vol;
        }

        Position pin_position( unsigned int ipin ) const {
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
            Position pos;
            pos.x = 0;
            pos.y = 0;
            pos.z = 0;

            int lat_x = ilat % nx_;
            int lat_y = ilat / nx_;

            for( int ix=0; ix<lat_x; ix++ ) {
                pos.x += this->at(ix, 0).nx();
            }
            pos.x += ipin % this->at(lat_x, lat_y).nx();

            for( int iy=0; iy<lat_y; iy++ ) {
                pos.y += this->at(0, iy).ny();
            }
            pos.y += ipin / this->at(lat_x, lat_y).ny();

            return pos;
        }

    private:
        // Plane dimensions in lattices
        unsigned int nx_;
        unsigned int ny_;

        unsigned int n_reg_;
        unsigned int n_xsreg_;
        
        // Locations of lattice interfaces
        VecF hx_;
        VecF hy_;
        // Local list of lattices
        std::vector<const Lattice*> lattices_;
        // List of the starting FSR index for each lattice in the plane
        VecI first_reg_lattice_;
    };
}
