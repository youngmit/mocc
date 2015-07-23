#pragma once

#define PI 3.1415926535897932
#define TWOPI 2.0*PI
#define HPI 0.5*PI
#define RPI 1.0/PI
#define RTWOPI 1.0/TWOPI
#define FPI 4.0*PI
#define RFPI 1.0/FPI


// Surface and direction indexing
enum Surface {
    EAST,
    NORTH,
    WEST,
    SOUTH,
    BOTTOM,
    TOP
};

// Boundary condition enumeration
enum Boundary {
    VACUUM,
    REFLECT,
    PARALLEL,
    PERIODIC
};

enum TraceDir {
    FW,
    BW
};
