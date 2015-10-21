#pragma once
#include <algorithm>
#include <vector>

#include "global_config.hpp"
#include "scattering_matrix.hpp"

namespace mocc{
    class Material{
    public:
    	Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                std::vector<VecF> scat);
        const VecF& xsab() const {
            return xsab_;
        }

        const VecF& xstr() const {
            return xstr_;
        }

        const VecF& xsnf() const {
            return xsnf_;
        }

        const VecF& xskf() const {
            return xsf_;
        }

        const VecF& xsch() const {
            return xsch_;
        }

        const ScatteringMatrix& xssc() const {
            return xssc_;
        }

        /**
         * Return whether the material is fissile.
         */
        bool is_fissile() const {
            return std::any_of(xsnf_.begin(), xsnf_.end(), 
                    [](real_t v){return v>0.0;});
        }

    private:
    	VecF xsab_;
    	VecF xstr_;
    	VecF xsnf_;
    	VecF xsf_;
    	VecF xsch_;
        ScatteringMatrix xssc_;
    };
}
