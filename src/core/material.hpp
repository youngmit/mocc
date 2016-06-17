/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once
#include <algorithm>
#include <vector>

#include "blitz_typedefs.hpp"
#include "global_config.hpp"
#include "scattering_matrix.hpp"

namespace mocc{
    class Material{
    public:
        Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch,
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

        // Accessors for fission
        const auto xsf() const {
            return xsf_;
        }
        const auto& xsf( int i ) const {
            return xsf_(i);
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
        ArrayB1 xsf_;
        ArrayB1 xsch_;
        ScatteringMatrix xssc_;
        bool is_fissile_;
    };
}
