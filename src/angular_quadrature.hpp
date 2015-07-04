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

        std::vector<Angle>::const_iterator begin() const {
            return angles_.cbegin();
        }

        std::vector<Angle>:: const_iterator end() const {
            return angles_.cend();
        }

        std::vector<Angle>::const_iterator octant( int octant ) const {
            assert( 0 < octant & octant < 9 );
            return angles_.cbegin() + (octant-1)*ndir_oct_;
        }

        ~AngularQuadrature() {
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
