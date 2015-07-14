#pragma once

#include <memory>
#include <vector>

#include "global_config.hpp"
#include "xs_mesh.hpp"
#include "source.hpp"

namespace mocc{
    class TransportSweeper{
    public:
        virtual ~TransportSweeper(){ }
        virtual void sweep(int group) =0;
        
        unsigned int n_reg() const {
            return n_reg_;
        }

        // Return a reference to the sweeper's XSMesh
        const XSMesh& xs_mesh() const {
            return xs_mesh_;
        }

        // Assign a source object to the sweeper.

    protected:
        unsigned int n_reg_;
        unsigned int ng_;

        XSMesh xs_mesh_;

        const Source* source_;

        VecF phis_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
