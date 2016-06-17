#pragma once

#include <iostream>

#define PI 3.1415926535897932
#define TWOPI (2.0*PI)
#define HPI (0.5*PI)
#define RPI (1.0/PI)
#define RTWOPI (1.0/TWOPI)
#define FPI (4.0*PI)
#define RFPI (1.0/FPI)

namespace mocc {
    // Surface and direction indexing
    enum Surface {
        EAST  = 0,
        NORTH = 1,
        WEST  = 2,
        SOUTH = 3,
        TOP = 4,
        BOTTOM = 5,
        INVALID = 6
    };



    enum class Direction : unsigned char {
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
         * Dirichlet boundary.
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
