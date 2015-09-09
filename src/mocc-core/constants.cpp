#include "constants.hpp"

namespace mocc {
extern const Surface AllSurfaces[] = { Surface::EAST, 
                                           Surface::NORTH, 
                                           Surface::WEST, 
                                           Surface::SOUTH, 
                                           Surface::TOP, 
                                           Surface::BOTTOM };

extern const Normal AllNormals[] = { Normal::X_NORM, Normal::Y_NORM,
    Normal::Z_NORM };
}
