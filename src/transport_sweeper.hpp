#pragma once

#include <memory>

namespace mocc{
    class TransportSweeper{
    public:
        virtual void sweep(int group) =0;
        
        int n_reg(){
            return n_reg_;
        }
    protected:
        int n_reg_;
    };

    typedef std::unique_ptr<TransportSweeper> UP_Sweeper_t;
}
