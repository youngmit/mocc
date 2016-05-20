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

#include "points.hpp"
#include "global_config.hpp"

namespace mocc {
    struct Particle {
    public:
        Particle() { }

        Particle(Point3 loc, Direction dir, int ig):
            location(loc),
            direction(dir),
            group(ig),
            weight(1.0)
        {
            return;
        }

        // Particle's location in the global domain
        Point3 location;
        Direction direction;

        // Particle's energy group
        int group;

        // Particle weight (for variance reduction and the like)
        real_t weight;
    };

} // namespace mocc
