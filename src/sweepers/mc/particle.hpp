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
#include <iosfwd>

#include "core/geometry/direction.hpp"
#include "core/geometry/points.hpp"
#include "core/global_config.hpp"
#include "position.hpp"

namespace mocc {
    /**
     * \brief Struct representing the state of a particle for Monte Carlo
     * simulation.
     *
     * The \ref Particle struct contains most of the state necessary to
     * represent a particle for Monte Carlo-simulation purposes. Right now, that
     * includes the particle's location in pin-local and global coordinates, its
     * direction of travel, and its energy group.
     *
     * There are a couple of hacks there, that should be cleaned up if the MC
     * stuff is to really go anywhere. First, the location is represented twice,
     * for global and pin-local coordinate, and the constructor only sets the
     * global coordinates. Moving to a more general, universe-based approach
     * would require more generality. Second, the location passed into the
     * constructor is assumed to be in global coordinates, and since we dont
     * want the Particle to have to be aware of the actual geometry it's moving
     * through, we need to manually set the pin-local position after
     * construction. Kind of weird, not too hard to fix, but good enough for
     * now. Another weird aspect of \ref Particle is that the global position is
     * 3-D while the local position is in 2-D. This is an artifact of the
     * underlying 2-D nature of the pin meshes. Again, to go more general would
     * require a change to the \ref Particle struct.
     */
    struct Particle {
    public:
        Particle() { }

        Particle(Point3 loc, Direction dir, int ig, unsigned id):
            weight(1.0),
            group(ig),
            direction(dir),
            location_global(loc),
            id(id),
            alive(true)
        {
            return;
        }
        // Particle weight (for variance reduction and the like)
        real_t weight;
        // Particle's energy group
        int group;
        // Particle's location in the pin-local domain
        Point2 location;
        Direction direction;
        // Particle's location in the global domain
        Point3 location_global;
        // Particles pin cell position
        Position pin_position;

        // ID, used for sorting and seeding RNG
        unsigned id;

        bool alive;

        /**
         * \brief Move the \ref Particle forward by a distance along its
         * direction of travel
         *
         * \param d the distance to move the \ref Particle
         *
         * This updates the \ref Particle's pin- and core-local positions
         */
        void move( real_t d ) {
            location.x += d*direction.ox;
            location.y += d*direction.oy;

            location_global.x += d*direction.ox;
            location_global.y += d*direction.oy;
            location_global.z += d*direction.oz;
            return;
        }

        /**
         * \brief A \ref Particle is considered "less than" another particle if
         * its ID is smaller.
         */
        bool operator<( const Particle &other ) const {
            return id<other.id;
        }

        friend std::ostream &operator<<(std::ostream &os, const Particle &p);
    };

} // namespace mocc
