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

#include <iosfwd>

namespace mocc {
struct Position {
    Position() : x(0), y(0), z(0)
    {
    }

    Position(unsigned int x, unsigned int y, unsigned int z) : x(x), y(y), z(z)
    {
    }

    int x;
    int y;
    int z;

    friend std::ostream &operator<<(std::ostream &os, const Position &pos);
};
}
