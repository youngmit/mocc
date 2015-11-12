#include "constants.hpp"

#include <error.hpp>

namespace mocc {
    const Surface AllSurfaces[] = { Surface::EAST,
                                    Surface::NORTH,
                                    Surface::WEST,
                                    Surface::SOUTH,
                                    Surface::TOP,
                                    Surface::BOTTOM };

    const Normal AllNormals[] = { Normal::X_NORM,
                                  Normal::Y_NORM,
                                  Normal::Z_NORM };

    std::ostream& operator<<(std::ostream& os, const Surface d ) {
        switch( d ) {
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
            case Surface::INVALID:
                os << "inv";
                break;
        }
        return os;
    }


    std::ostream& operator<<(std::ostream& os, const Direction d ) {
        switch( d ) {
            case Direction::EAST:
                os << "east";
                break;
            case Direction::WEST:
                os << "west";
                break;
            case Direction::NORTH:
                os << "north";
                break;
            case Direction::SOUTH:
                os << "south";
                break;
            case Direction::TOP:
                os << "top";
                break;
            case Direction::BOTTOM:
                os << "bottom";
                break;
            case Direction::NE:
                os << "ne";
                break;
            case Direction::NW:
                os << "nw";
                break;
            case Direction::SW:
                os << "sw";
                break;
            case Direction::SE:
                os << "se";
                break;
            case Direction::INVALID:
                os << "inv";
                break;
        }
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const Normal n ) {
        switch( n ) {
            case Normal::X_NORM:
                os << "X-Normal";
                break;
            case Normal::Y_NORM:
                os << "Y-Normal";
                break;
            case Normal::Z_NORM:
                os << "Z-Normal";
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
