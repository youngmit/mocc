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

    private:
        // Enumerated quadrature type
        quad_t type_;

        // Number of angles per octant
        int ndir_oct_;

        // Vector of all angles for all octants
        std::vector<Angle> angles_;

    };
}
