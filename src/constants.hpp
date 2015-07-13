#pragma once

#define PI 3.1415926535897932
#define TWOPI 2.0*PI
#define HPI 0.5*PI
#define RPI 1.0/PI
#define RTWOPI 1.0/TWOPI


// Surface and direction indexing
enum Surface {
    EAST,
    NORTH,
    WEST,
    SOUTH,
    BOTTOM,
    TOP
};

enum TraceDir {
    FW,
    BW
};
