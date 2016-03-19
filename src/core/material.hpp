#pragma once
#include <algorithm>
#include <vector>

#include "blitz_typedefs.hpp"
#include "global_config.hpp"
#include "scattering_matrix.hpp"

namespace mocc{
    class Material{
    public:
        Material(VecF xsab, VecF xsnf, VecF xskf, VecF xsch,
                std::vector<VecF> scat);
        // Accessors for absorption
        const auto xsab() const {
            return xsab_;
        }
        const auto& xsab( int i ) const {
            return xsab_(i);
        }

        // Accessors for transport xs
        const auto xstr() const {
            return xstr_;
        }
        const auto& xstr( int i ) const {
            return xstr_(i);
        }

        // Accessors for nu-fission
        const auto xsnf() const {
            return xsnf_;
        }
        const auto& xsnf( int i ) const {
            return xsnf_(i);
        }

        // Accessors for kappa-fission
        const auto xskf() const {
            return xskf_;
        }
        const auto& xskf( int i ) const {
            return xskf_(i);
        }

        // Accessors for chi
        const auto xsch() const {
            return xsch_;
        }
        const auto& xsch( int i ) const {
            return xsch_(i);
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
        ArrayB1 xsab_;
        ArrayB1 xstr_;
        ArrayB1 xsnf_;
        ArrayB1 xskf_;
        ArrayB1 xsch_;
        ScatteringMatrix xssc_;
        bool is_fissile_;
    };
}
