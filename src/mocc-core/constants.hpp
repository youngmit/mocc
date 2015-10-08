#pragma once

#include <iostream>

#define PI 3.1415926535897932
#define TWOPI 2.0*PI
#define HPI 0.5*PI
#define RPI 1.0/PI
#define RTWOPI 1.0/TWOPI
#define FPI 4.0*PI
#define RFPI 1.0/(FPI)

namespace mocc {
    // Surface and direction indexing
    enum Surface {
        EAST  = 0,
        NORTH = 1,
        WEST  = 2,
        SOUTH = 3,
        TOP = 4,
        BOTTOM = 5,
        INVALID = 6,
    };



    enum class Direction {
        EAST  = 0,
        NORTH = 1,
        WEST  = 2,
        SOUTH = 3,
        TOP = 4,
        BOTTOM = 5,
        NE,
        NW,
        SW,
        SE,
        INVALID
    };

    extern const Surface AllSurfaces[6];
        
    enum class Normal {
        X_NORM = 0,
        Y_NORM,
        Z_NORM
    };

    extern const Normal AllNormals[3];
    
    
    // Boundary condition enumeration
    enum class Boundary {
        VACUUM,
        REFLECT,
        PARALLEL,
        PERIODIC,
        INVALID
    };
    
    enum class TraceDir {
        FW,
        BW
    };

    std::ostream& operator<<(std::ostream& os, const Surface s );

    Normal surface_to_normal( Surface s );
    
}
