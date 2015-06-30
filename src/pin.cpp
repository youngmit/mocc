#include "pin.hpp"

namespace mocc {
    Pin::Pin( int id, PinMesh* pin, VecI mat ):
        id_(id),
        pin_mesh_(pin),
        mat_IDs_(mat) {}
}
