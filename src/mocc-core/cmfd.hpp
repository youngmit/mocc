#pragma once

#include <memory>

#include "global_config.hpp"
#include "coarse_data.hpp"
#include "eigen_interface.hpp"
#include "mesh.hpp"

namespace mocc {
    class CMFD {
    public:
        CMFD( Mesh *mesh, int ng );

        void solve();

        void project( ArrayX &flux );

        /** 
         * Return a pointer to the coarse data. This is used to couple sweepers
         * and other objects that need access to the coarse data to the CMFD
         * solver.
         */
        CoarseData* get_data(){
            return &coarse_data_;
        }
    private:
        Mesh* mesh_;
        CoarseData coarse_data_;
    };
    typedef std::unique_ptr<CMFD> UP_CMFD_t;
}
