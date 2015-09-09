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
}
