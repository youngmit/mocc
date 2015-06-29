#pragma once

#include <memory>

namespace mocc{
class TransportSweeper{
public:
    virtual void sweep() =0;
    
    int n_reg(){
        return n_reg_;
    }
protected:
    int n_reg_;
};

}
