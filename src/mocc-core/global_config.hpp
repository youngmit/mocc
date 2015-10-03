#pragma once

#include <cstddef>
#include <vector>
#include <valarray>
#include <iostream>

namespace mocc{

#ifdef FORCE_SINGLE
typedef float float_t;
#define FLOAT_EPS 1e-5
#else
typedef double float_t;
#define FLOAT_EPS 1e-12
#endif

// General purpose vector of floats, ints, etc
typedef std::vector<float_t> VecF;
typedef std::vector<unsigned int> VecI;
typedef std::vector<int> VecSI;

// valarray of floats
typedef std::valarray<float_t> ArrayF;


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
    unsigned int x;
    unsigned int y;
    unsigned int z;

    friend std::ostream& operator<<(std::ostream& os, const Position &pos ) {
        os << pos.x << " " << pos.y << " " << pos.z;
        return os;
    }
};

}

#define PROG_NAME "MOCC"

