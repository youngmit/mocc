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

#include "fission_bank.hpp"

#include "pugixml.hpp"

#include "error.hpp"

namespace mocc {
    FissionBank::FissionBank():
        mesh_(nullptr)
    {
        return;
    }

    /**
     * \brief Construct a FissionBank by uniformly sampling fission sites.
     *
     * \param input XML node containing bounds of a 3-D box within which to
     * sample initial fission sites
     * \param n the number of initial sites to sample
     *
     * This constructor initializes a \ref FissionBank using input specified in
     * an XML node.
     */
    FissionBank::FissionBank( const pugi::xml_node &input, int n,
            const CoreMesh &mesh ):
        mesh_(&mesh)
    {
        if( input.empty() ) {
            throw EXCEPT("Empty input provided to FissionBank");
        }

        // Make sure that all of the bounds are specified.
        if( input.attribute("x_min").empty() ||
            input.attribute("x_max").empty() ||
            input.attribute("y_min").empty() ||
            input.attribute("y_max").empty() ||
            input.attribute("z_min").empty() ||
            input.attribute("z_max").empty() )
        {
            throw EXCEPT("Not all X, Y, Z bounds specified in fission_box");
        }

        real_t x_min = input.attribute("x_min").as_double(0.0);
        real_t x_max = input.attribute("x_max").as_double(0.0);
        real_t y_min = input.attribute("y_min").as_double(0.0);
        real_t y_max = input.attribute("y_max").as_double(0.0);
        real_t z_min = input.attribute("z_min").as_double(0.0);
        real_t z_max = input.attribute("z_max").as_double(0.0);

        // Make sure the bounds are valid
        if( (x_min >= x_max) || (y_min >= y_max) || (z_min >= z_max) ) {
            throw EXCEPT("Invalid fission_box bounds specified.");
        }

        // See if we want to do a fissile region rejection (only accept fission
        // sites in fissile regions).
        bool fissile_rejection = input.attribute("fissile_rejection").
            as_bool(true);

        sites_.reserve( n );
        if( !fissile_rejection ) {
            for( int i=0; i<n; i++ ) {
                sites_.push_back({ RNG_MC.random(x_min, x_max),
                                   RNG_MC.random(y_min, y_max),
                                   RNG_MC.random(z_min, z_max) });
            }
        } else {
            Warn("Fissile region rejection is not supported yet.");
            for( int i=0; i<n; i++ ) {
                sites_.push_back({ RNG_MC.random(x_min, x_max),
                                   RNG_MC.random(y_min, y_max),
                                   RNG_MC.random(z_min, z_max) });
            }
        }

        return;
    }

    real_t FissionBank::shannon_entropy() const {
        real_t h = 0.0;

        throw EXCEPT("not implemented");

        return h;
    }

    void FissionBank::swap( FissionBank &other ) {
        sites_.swap( other.sites_ );
    }
} // namespace mocc
