#include "material.hpp"

#include <algorithm>

namespace mocc{
    Material::Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch,
                       std::vector<VecF> scat):
        xssc_(scat)
    {
        xsab_ = xsab;
        xsnf_ = xsnf;
        xsf_  = xsf;
        xsch_ = xsch;

        // Normalize chi. Wow. Very precision. Much digits.
        if( std::any_of(xsch_.begin(), xsch_.end(),
                    [](real_t v){return v>0.0;}) ) {
            real_t chi_sum = 0.0;
            for( auto c: xsch_ ) {
                chi_sum += c;
            }
            for( auto &c: xsch_ ) {
                c /= chi_sum;
            }
        }

        int ng = xsab_.size();
        xstr_ = VecF( ng, 0.0 );
        // Simple calculation of transport cross section
        for( int ig=0; ig<ng; ig++ ) {
            xstr_[ig] = xsab[ig] + xssc_.out(ig);
        }
    }
};

