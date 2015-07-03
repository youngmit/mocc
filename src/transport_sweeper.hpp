#pragma once

#include <memory>
#include <vector>

#include "global_config.hpp"

namespace mocc{
    class TransportSweeper{
    public:
        virtual ~TransportSweeper(){ }
        virtual void sweep(int group) =0;
        
        int n_reg(){
            return n_reg_;
        }
    protected:
        int n_reg_;
        int ng_;

        VecF phis_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
