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

#include "constants.hpp"

#include <iostream>
#include "util/error.hpp"

namespace mocc {
const Surface AllSurfaces[] = {Surface::EAST,  Surface::NORTH, Surface::WEST,
                               Surface::SOUTH, Surface::TOP,   Surface::BOTTOM};

const Normal AllNormals[] = {Normal::X_NORM, Normal::Y_NORM, Normal::Z_NORM};

std::ostream &operator<<(std::ostream &os, const Surface d)
{
    switch (d) {
    case Surface::EAST:
        os << "east";
        break;
    case Surface::WEST:
        os << "west";
        break;
    case Surface::NORTH:
        os << "north";
        break;
    case Surface::SOUTH:
        os << "south";
        break;
    case Surface::TOP:
        os << "top";
        break;
    case Surface::BOTTOM:
        os << "bottom";
        break;
    case Surface::INTERNAL:
        os << "internal";
        break;
    case Surface::INVALID:
        os << "inv";
        break;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Cardinal d)
{
    switch (d) {
    case Cardinal::EAST:
        os << "east";
        break;
    case Cardinal::WEST:
        os << "west";
        break;
    case Cardinal::NORTH:
        os << "north";
        break;
    case Cardinal::SOUTH:
        os << "south";
        break;
    case Cardinal::TOP:
        os << "top";
        break;
    case Cardinal::BOTTOM:
        os << "bottom";
        break;
    case Cardinal::NE:
        os << "ne";
        break;
    case Cardinal::NW:
        os << "nw";
        break;
    case Cardinal::SW:
        os << "sw";
        break;
    case Cardinal::SE:
        os << "se";
        break;
    case Cardinal::INVALID:
        os << "inv";
        break;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Normal n)
{
    switch (n) {
    case Normal::X_NORM:
        os << "X-Normal";
        break;
    case Normal::Y_NORM:
        os << "Y-Normal";
        break;
    case Normal::Z_NORM:
        os << "Z-Normal";
        break;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Boundary b)
{
    switch (b) {
    case Boundary::REFLECT:
        os << "REFLECT";
        break;
    case Boundary::VACUUM:
        os << "VACUUM";
        break;
    case Boundary::PARALLEL:
        os << "PARALLEL";
        break;
    case Boundary::PERIODIC:
        os << "PERIODIC";
        break;
    case Boundary::PRESCRIBED:
        os << "PRESCRIBED";
        break;
    default:
        os << "Unknown: " << (int)b;
    }

    return os;
}

Normal surface_to_normal(Surface s)
{
    switch (s) {
    case Surface::EAST:
        return Normal::X_NORM;
    case Surface::WEST:
        return Normal::X_NORM;
    case Surface::NORTH:
        return Normal::Y_NORM;
    case Surface::SOUTH:
        return Normal::Y_NORM;
    case Surface::BOTTOM:
        return Normal::Z_NORM;
    case Surface::TOP:
        return Normal::Z_NORM;
    default:
        throw EXCEPT("Unsupported surface.");
    }
}
}
