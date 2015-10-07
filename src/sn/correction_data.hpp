#pragma once

#include <memory>

#include "constants.hpp"
#include "global_config.hpp"

namespace mocc {

    /**
     * This class provides a storage scheme for the correction factors needed to
     * perform corrected diamond difference. The CDD Sn and MoC sweepers must be
     * provided with a reference to an object of this class to access and store
     * correction factors, respectively. Due to the relatively high
     * dimensionality of the data (space, angle, energy and cardinal direction
     * [X|Y]), instead of using a multidimensional array, we will instead use
     * accessor functions to get the data out of a dense linear representation.
     *
     * The factors are stored 
     */
    class CorrectionData {
    public:
        CorrectionData( ):
            nreg_( 0 ),
            nang_( 0 ),
            ngroup_( 0 )
        {
            return;
        }

        CorrectionData( size_t nreg, size_t nang, size_t ngroup ):
            nreg_( nreg ),
            nang_( nang ),
            ngroup_( ngroup ),
            alpha_( nreg_*nang_*ngroup_*2, 0.5 ),
            beta_( nreg_*nang_*ngroup_, 1.0 )
        {
            return;
        }
            

        real_t& alpha( size_t reg, size_t ang, size_t group, Normal norm ) {
            return alpha_[ nreg_*nang_*2*group + nreg_*2*ang + 2*reg + 
                (int)norm ];
        }
        const real_t& alpha( size_t reg, size_t ang, size_t group, 
                Normal norm ) const {
            return alpha_[ nreg_*nang_*2*group + nreg_*2*ang + 2*reg + 
                (int)norm ];
        }

        real_t& beta( size_t reg, size_t ang, size_t group ) {
            return beta_[ nreg_*nang_*group + nreg_*ang + reg ];
        }

        const real_t& beta( size_t reg, size_t ang, size_t group ) const {
            return beta_[ nreg_*nang_*group + nreg_*ang + reg ];
        }

    private:
        size_t nreg_;
        size_t nang_;
        size_t ngroup_;

        VecF alpha_;
        VecF beta_;
    };

    typedef std::unique_ptr<CorrectionData> UP_CorrectionData_t;
}
