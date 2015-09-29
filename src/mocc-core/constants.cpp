#include "constants.hpp"

#include <error.hpp>

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
        os << "(" << (int)s << ") ";
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

    Normal surface_to_normal( Surface s ) {
        switch( s ) {
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
                std::cout << (int)s << std::endl;
                throw EXCEPT("Unsupported surface.");
        }
    }
}
