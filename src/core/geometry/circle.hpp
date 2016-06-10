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

#include "core/global_config.hpp"

#include "direction.hpp"
#include "points.hpp"

namespace mocc {
struct Circle {
    Circle(Point2 c, real_t r) : c(c), r(r)
    {
    }

    Point2 c;
    real_t r;

    /// \todo document this
    real_t distance_to_surface(Point2 p, Direction dir) const;
};

} // namespace mocc
