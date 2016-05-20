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

#pragma once

#include <array>
#include <cassert>
#include <map>
#include <vector>

#include "assembly.hpp"
#include "constants.hpp"
#include "geometry/geom.hpp"
#include "pugifwd.hpp"

namespace mocc {
    class Core {
    public:
        Core();
        Core( const pugi::xml_node &input,
              const std::map<int, UP_Assembly_t> &assemblies);
        ~Core();

        const Assembly& at( unsigned int i ) const {
            assert( (0 <= i) & (i < assemblies_.size()) );
            return *assemblies_[i];
        }

        const Assembly& at( unsigned int x, unsigned int y ) const {
            assert( (0 <= x) & (x < nx_) );
            assert( (0 <= y) & (y < ny_) );
            return *assemblies_[y*nx_ + x];
        }

        const std::vector<Assembly*>& assemblies() const {
            return assemblies_;
        }

        /**
         * Return the number of assemblies along X direction
         */
        int nx() const {
            return nx_;
        }

        /**
         * Return the number of assemblies along Y direction
         */
        int ny() const {
            return ny_;
        }

        /**
         * Return the total number of assemblies in the core
         */
        int nasy() const {
            return assemblies_.size();
        }

        /**
         * Return the number of pins along the X direction
         */
        int npin_x() const {
            return npinx_;
        }

        /**
         * Return the number of pins along the Y direction
         */
        int npin_y() const {
            return npiny_;
        }

        /**
         * Return the number of planes in the core
         */
        int nz() const {
            return assemblies_[0]->nz();
        }

        /**
         * Return the plane heights
         */
        const VecF& dz() const {
            return assemblies_.front()->dz();
        }

        /**
         * Return the boundary condition array
         */
        std::array<Boundary, 6> boundary() const {
            return bc_;
        }

    private:
        // Core dimensions (in assemblies)
        unsigned int nx_;
        unsigned int ny_;

        // Core dimensions (in pins)
        unsigned int npinx_;
        unsigned int npiny_;

        // Boundaries of the lattices in the core
        VecF hx_vec_;
        VecF hy_vec_;


        // 2D array of assemblies
        std::vector<Assembly*> assemblies_;

        // boundary conditions
        std::array<Boundary, 6> bc_;
    };

    Core ParseCore( const pugi::xml_node &input,
            const std::map<int, UP_Assembly_t> &assemblies );
}
