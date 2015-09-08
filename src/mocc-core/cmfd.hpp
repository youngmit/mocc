#pragma once

#include "global_config.hpp"
#include "eigen_interface.hpp"
#include "mesh.hpp"

namespace mocc {
    class CMFD {
    public:
        CMFD( ) {

        }

        CMFD( Mesh *mesh );

        void solve();

        void project( ArrayX &flux );
    private:
        Mesh* mesh_;
    };
}
