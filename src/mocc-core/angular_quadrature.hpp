#pragma once

#include <cassert>
#include <iostream>
#include <vector>

#include "pugixml.hpp"

#include "angle.hpp"
#include "global_config.hpp"

namespace mocc {

    enum quad_t {
        SN // Level-symmetric
    };


    class AngularQuadrature {
    public:
        AngularQuadrature( const pugi::xml_node &input );

        ~AngularQuadrature() {
        }

        // Return an iterator to the first angle in the quadrature
        std::vector<Angle>::const_iterator begin() const {
            return angles_.cbegin();
        }

        // Return an iterator past the last angle in the quadrature
        std::vector<Angle>:: const_iterator end() const {
            return angles_.cend();
        }

        /**
         * 
         * Return an iterator to the first angle in the given octant. Octants
         * are indexed from 1, following mathematical convention. Also,
         * following convention for container classes, specifying octant 9, is
         * tantamount to end().
         */
        std::vector<Angle>::const_iterator octant( int octant ) const {
            return angles_.cbegin() + (octant-1)*ndir_oct_;
        }

        /*
         * Return a const reference to the angle indexed
         */
        const Angle& operator[]( int iang ) {
            return angles_[iang];
        }

        /**
         * Return the number of angles in each octant
         */
        int ndir_oct() const {
            return ndir_oct_;
        }

        /**
         * Return the total number of angles
         */
        int ndir() const {
            return angles_.size();
        }

        /**
         * Modify one of the angles in the quadrature. The new angle provided
         * should be specified on the first octant, and all corresponding angles
         * in other octants are updated internally.
         */
        void modify_angle( int iang, Angle ang );

        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, 
                const AngularQuadrature &angquad) {
            os << "Alpha\tTheta\tomega x   \tomega y   \tomega z" << std::endl;
            for( auto &ang: angquad.angles_ ) {
                os << ang << std::endl;
            }
            return os;
        }

        /**
         * Return the index of the angle reflected across a surface with the
         * given normal.
         */
        unsigned int reflect( unsigned int iang, Normal normal ) const {
            int ioct = iang / ndir_oct_;
            int new_oct = 0;
            // Lets be real, im just showing off here...
            switch(normal) {
                case Normal::X_NORM:
                    new_oct = ioct + 1-2*(ioct%2);
                    break;
                case Normal::Y_NORM:
                    new_oct = abs(ioct + 3-(ioct%2)*6) % 4;
                    break;
                case Normal::Z_NORM:
                    new_oct = (ioct+4) % 8;
                    break;
            }

            return iang + (new_oct - ioct)*ndir_oct_;
        }

        /**
         * Return the index of the angle reflected across the given surface
         */
        unsigned int reflect( unsigned int iang, Surface surf ) const {
            if( (surf == Surface::NORTH) | (surf == Surface::SOUTH) ) {
                return this->reflect( iang, Normal::Y_NORM );
            } else if ( (surf == Surface::EAST) | (surf == Surface::WEST) ) {
                return this->reflect( iang, Normal::X_NORM );
            } else {
                return this->reflect( iang, Normal::Z_NORM );
            }

        }

        /**
         * Return the index of the angle that is in the reverse direction of the
         * angle index passed. This can operate in two different modes, based on
         * dim, which should be 2[D] or 3[D]. For 2D, the returned angle always
         * lies in the positive-Z half-space. For 3D, the returned angle
         */
        unsigned int reverse( unsigned int iang, unsigned int dim=2 ) const {
            assert( (dim == 2) || (dim ==3) );
            if( dim == 2) {
                return (iang + ndir_oct_*2) % (ndir_oct_*4);
            } else {
                return (iang + ndir_oct_*6) % (ndir_oct_*8);
            }
            return 0;
        }
    private:
        // Enumerated quadrature type
        quad_t type_;

        // Number of angles per octant
        int ndir_oct_;

        // Vector of all angles for all octants
        std::vector<Angle> angles_;

    };
}
