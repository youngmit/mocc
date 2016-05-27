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

#define PI 3.1415926535897932
#define TWOPI (2.0*PI)
#define HPI (0.5*PI)
#define RPI (1.0/PI)
#define RTWOPI (1.0/TWOPI)
#define FPI (4.0*PI)
#define RFPI (1.0/FPI)

namespace mocc {
    // Surface and direction indexing
    // This is not an enum class, because it is used elsewhere in a bitfield (i
    // know, i know), and compilers tend to complain when you try to guarantee
    // to them in one place that they get 8 bits and another that they get 4 or
    // whatever.
    enum Surface {
        EAST  = 0,
        NORTH = 1,
        WEST  = 2,
        SOUTH = 3,
        TOP = 4,
        BOTTOM = 5,
        INTERNAL = 6,
        INVALID = 7
    };

    enum class Reaction : unsigned char {
        SCATTER = 0,
        FISSION = 1,
        CAPTURE = 2
    };

    enum class Cardinal : unsigned char {
        EAST  = 0,
        NORTH = 1,
        WEST  = 2,
        SOUTH = 3,
        TOP = 4,
        BOTTOM = 5,
        NE = 6,
        NW = 7,
        SW = 8,
        SE = 9,
        INVALID = 10
    };

    extern const Surface AllSurfaces[6];

    enum class Normal : unsigned char {
        X_NORM = 0,
        Y_NORM = 1,
        Z_NORM = 2
    };

    extern const Normal AllNormals[3];

    // Boundary condition enumeration
    enum class Boundary {
        /**
         * Zero incoming flux
         */
        VACUUM,
        /**
         * Reflected incoming flux
         */
        REFLECT,
        /**
         * Incoming flux communicated between domain nodes
         */
        PARALLEL,
        /**
         * Flux exiting one face enters the opposite face, same angle
         */
        PERIODIC,
        /**
         * Boundary condition prescribed as incoming angular flux, namely,
         * Dirichelet boundary.
         */
        PRESCRIBED,

        INVALID
    };

    enum class TraceDir {
        FW,
        BW
    };

    std::ostream& operator<<(std::ostream& os, const Surface s );

    std::ostream& operator<<(std::ostream& os, const Boundary b );

    std::ostream& operator<<(std::ostream& os, const Normal n );

    Normal surface_to_normal( Surface s );

}
