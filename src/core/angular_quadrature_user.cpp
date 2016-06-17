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
#include "angular_quadrature_user.hpp"

#include "pugixml.hpp"

#include "error.hpp"

namespace mocc {

    std::vector<Angle> GenUserQuadrature( const pugi::xml_node &input ) {
        std::vector<Angle> angles;

        // Count the number of specified angles
        int n = 0;
        for( auto node = input.child("angle"); node;
                node = node.next_sibling("angle") ) {
            n++;
        }

        angles.reserve(n);

        // Loop back through the angles and read them in. For the most part, we
        // will just rely on the Angle constructor that uses XML input. Slick!
        for( auto node = input.child("angle"); node;
                node = node.next_sibling("angle") ) {
            angles.push_back( Angle(node) );
        }

        // make sure all of the angles are in octant 1, and that their weights
        // sum to 1. The AngularQuadrature constructor is going to try and
        // expand them to all angles, so this is important.
        real_t w = 0.0;
        for( const auto &angle: angles ) {
            w += angle.weight;
            if( (angle.ox < 0.0) || (angle.oy < 0.0) || angle.oz < 0.0 ) {
                throw EXCEPT("User-specified angle is not in octant 1.");
            }
        }

        // This might need to be relaxed to not be super annoying. Perhaps allow
        // more variation from unity, but scale the weights to unity within
        // machine precision.
        if( !fp_equiv_ulp(w, 1.0) ) {
            throw EXCEPT("User-specified angle weights do not sum to one in "
                    "first octant");
        }

        return angles;
    }

} // namespace mocc
