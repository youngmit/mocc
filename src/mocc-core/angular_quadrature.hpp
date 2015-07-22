#pragma once

#include <vector>
#include <cassert>
#include <iostream>

#include "pugixml.hpp"

#include "global_config.hpp"
#include "angle.hpp"

namespace mocc {

    enum quad_t {
        SN // Level-symmetric
    };


    class AngularQuadrature {
    public:
        AngularQuadrature( const pugi::xml_node &input );

        ~AngularQuadrature() {
        }

        std::vector<Angle>::const_iterator begin() const {
            return angles_.cbegin();
        }

        std::vector<Angle>:: const_iterator end() const {
            return angles_.cend();
        }

        std::vector<Angle>::const_iterator octant( int octant ) const {
            return angles_.cbegin() + (octant-1)*ndir_oct_;
        }

        // Return the number of angles in each octant
        int ndir_oct() const {
            return ndir_oct_;
        }

        // Modify one of the angles in the quadrature. The new angle provided
        // should be specified on the first octant, and all corresponding angles
        // in other octants are updated internally,
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

        // Return the index of the angle reflected across the given surface
        unsigned int reflect( unsigned int iang, Surface surf ) const {
            int ioct = iang / ndir_oct_;
            assert( ioct < 4 ); // Not doing negative-z angles for now
            int new_oct = 0;
            // Lets be real, im just showing off here...
            if( (surf == NORTH) | (surf == SOUTH) ) {
                new_oct = abs(ioct + 3-(ioct%2)*6) % 4;
            } else if ( (surf == EAST) | (surf == WEST) ) {
                new_oct = ioct + 1-2*(ioct%2);
            }

            return iang + (new_oct - ioct)*ndir_oct_;
        }

        // Return the index of the angle that is in the reverse direction of the
        // angle index passed. This can operate in two different modes, based on
        // dim, which should be 2[D] or 3[D]. For 2D, the returned angle always
        // lies in the positive-Z half-space. For 3D, the returned angle 
        unsigned int reverse( unsigned int iang, unsigned int dim=2 ) const {
            assert( dim == 2 | dim ==3 );
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
