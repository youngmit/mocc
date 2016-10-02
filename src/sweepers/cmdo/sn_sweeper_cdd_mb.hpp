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

#include "sn_sweeper_cdd.hpp"

namespace mocc {
namespace cmdo {
/**
 * \brief Specialization of the CDD sweeper for Primitive Multiple-Balance
 */
class SnSweeper_CDD_PMB : public SnSweeper_CDD<SnSweeper_CDD_PMB> {
public:
    SnSweeper_CDD_PMB(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeper_CDD<SnSweeper_CDD_PMB>(input, mesh)
    {
        return;
    }

    real_t evaluate_2d(real_t &flux_x, real_t &flux_y, real_t q, real_t xstr,
                       int i, const ThreadState &t_state) const
    {
        return SnSweeper_CDD<SnSweeper_CDD_PMB>::evaluate_2d(flux_x, flux_y, q,
                                                             xstr, i, t_state);
    }

    real_t evaluate(real_t &flux_x, real_t &flux_y, real_t &flux_z, real_t q,
                    real_t xstr, int i, const ThreadState &t_state) const
    {
        int ix = i % mesh_.nx();
        int ia = t_state.macroplane * this->plane_size_ + i % this->plane_size_;
        real_t tx = t_state.ox / mesh_.dx(ix);

        real_t ax =
            corrections_->alpha(ia, t_state.iang_2d, group_, Normal::X_NORM);
        real_t ay =
            corrections_->alpha(ia, t_state.iang_2d, group_, Normal::Y_NORM);
        real_t b = corrections_->beta(ia, t_state.iang_2d, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t psi = q - (q / (2 + (xstr / t_state.tz))) +
                     2.0 * (tx * flux_x + t_state.ty * flux_y) +
                     t_state.tz * flux_z;
        psi /=
            tx / gx + t_state.ty / gy + t_state.tz / (1.0 + xstr * 0.5) + xstr;

        flux_x = psi / gx - flux_x;
        flux_y = psi / gy - flux_y;
        flux_z = (2.0 * psi * t_state.tz + q) / (2.0 * t_state.tz + xstr);

        return psi;
    }
};
}
}
