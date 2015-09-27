#include "constants.hpp"

namespace mocc {
const Surface AllSurfaces[] = { Surface::EAST, 
                                           Surface::NORTH, 
                                           Surface::WEST, 
                                           Surface::SOUTH, 
                                           Surface::TOP, 
                                           Surface::BOTTOM };

const Normal AllNormals[] = { Normal::X_NORM, Normal::Y_NORM,
    Normal::Z_NORM };


std::ostream& operator<<(std::ostream& os, const Surface s ) {
    switch( s ) {
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
        case Surface::NE:
            os << "ne";
            break;
        case Surface::NW:
            os << "nw";
            break;
        case Surface::SW:
            os << "sw";
            break;
        case Surface::SE:
            os << "se";
            break;
        case Surface::INVALID:
            os << "invalid";
            break;
    }
    return os;
}
}
