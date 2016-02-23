#pragma once

#include <cstddef>
#include <vector>
#include <valarray>
#include <iostream>

namespace mocc{

#ifdef FORCE_SINGLE
typedef float real_t;
#define FLOAT_EPS 1e-5
#else
typedef double real_t;
#define FLOAT_EPS 1e-12
#endif

// General purpose vector of floats, ints, etc
typedef std::vector<real_t> VecF;
typedef std::vector<int> VecI;
typedef std::vector<int> VecSI;

// valarray of floats
typedef std::valarray<real_t> ArrayF;
typedef std::slice_array<real_t> SliceF;


struct Position {
    Position():
        x( 0 ),
        y( 0 ),
        z( 0 )
    { }

    Position( unsigned int x, unsigned int y, unsigned int z):
        x( x ),
        y( y ),
        z( z )
    { }
    int x;
    int y;
    int z;

    friend std::ostream& operator<<(std::ostream& os, const Position &pos ) {
        os << pos.x << " " << pos.y << " " << pos.z;
        return os;
    }
};

}

#define PROG_NAME "MOCC"

