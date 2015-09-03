#pragma once

#define PI 3.1415926535897932
#define TWOPI 2.0*PI
#define HPI 0.5*PI
#define RPI 1.0/PI
#define RTWOPI 1.0/TWOPI
#define FPI 4.0*PI
#define RFPI 1.0/(FPI)


// Surface and direction indexing
enum class Surface {
    EAST = 0,
    NORTH,
    WEST,
    SOUTH,
    TOP,
    BOTTOM,
    NE,
    NW,
    SW,
    SE,
    INVALID
};
static const Surface AllSurfaces[] = { Surface::EAST, 
                                       Surface::NORTH, 
                                       Surface::WEST, 
                                       Surface::SOUTH, 
                                       Surface::TOP, 
                                       Surface::BOTTOM };

enum class Normal {
    X_NORM = 0,
    Y_NORM,
    Z_NORM
};
static Normal AllNormals[] = { Normal::X_NORM, Normal::Y_NORM, Normal::Z_NORM };

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
