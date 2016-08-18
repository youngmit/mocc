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

#include <iostream>

#include "util/force_inline.hpp"
#include "core/mesh.hpp"
#include "sn/cell_worker.hpp"
#include "sn/sn_sweeper.hpp"

namespace mocc {
namespace sn {
/**
 * A simple class extending \ref sn::CellWorker to perform the algebraic
 * work needed to propagate flux through an orthogonal mesh cell using
 * the diamond difference scheme.
 */
class SnSweeper_DD : public SnSweeperVariant<SnSweeper_DD> {
public:
    SnSweeper_DD(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeperVariant<SnSweeper_DD>(input, mesh)
    {
        return;
    }

    real_t MOCC_FORCE_INLINE evaluate(real_t &flux_x, real_t &flux_y,
                                      real_t &flux_z, real_t q, real_t xstr,
                                      int i, const ThreadState &t_state) const
    {
        size_t ix = i % mesh_.nx();
        real_t tx = t_state.angle.ox / mesh_.dx(ix);
        real_t psi =
            2.0 * (tx * flux_x + t_state.ty * flux_y + t_state.tz * flux_z) + q;
        psi /= 2.0 * (tx + t_state.ty + t_state.tz) + xstr;

        flux_x = 2.0 * psi - flux_x;
        flux_y = 2.0 * psi - flux_y;
        flux_z = 2.0 * psi - flux_z;


        return psi;
    }

    real_t evaluate_2d(real_t &flux_x, real_t &flux_y, real_t q, real_t xstr,
                       int i, const ThreadState &t_state) const
    {
        size_t ix  = i % mesh_.nx();
        real_t tx  = t_state.angle.ox / mesh_.dx(ix);
        real_t psi = 2.0 * (tx * flux_x + t_state.ty * flux_y) + q;
        psi /= 2.0 * (tx + t_state.ty) + xstr;

        flux_x = 2.0 * psi - flux_x;
        flux_y = 2.0 * psi - flux_y;

        return psi;
    }
};

class SnSweeper_DD_SC : public SnSweeperVariant<SnSweeper_DD_SC> {
public:
    SnSweeper_DD_SC(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeperVariant<SnSweeper_DD_SC>(input, mesh)
    {
        return;
    }

    real_t evaluate_2d(real_t &flux_x, real_t &flux_y, real_t q, real_t xstr,
                       int i, const ThreadState &t_state) const
    {
        size_t ix  = i % mesh_.nx();
        real_t tx  = t_state.angle.ox / mesh_.dx(ix);
        real_t psi = 2.0 * (tx * flux_x + t_state.ty * flux_y) + q;
        psi /= 2.0 * (tx + t_state.ty) + xstr;

        flux_x = 2.0 * psi - flux_x;
        flux_y = 2.0 * psi - flux_y;

        return psi;
    }

    real_t evaluate(real_t &flux_x, real_t &flux_y, real_t &flux_z, real_t q,
                    real_t xstr, int i, const ThreadState &t_state) const
    {
        size_t ix = i % mesh_.nx();
        real_t tx = t_state.angle.ox / mesh_.dx(ix);

        real_t tau    = xstr / t_state.tz;
        real_t rho    = 1.0 / tau - 1.0 / (std::exp(tau) - 1.0);
        real_t rhofac = rho / (1.0 - rho);

        real_t psi = 2.0 * (tx * flux_x + t_state.ty * flux_y) +
                     t_state.tz * (rhofac + 1.0) * flux_z + q;
        psi /= 2.0 * (tx + t_state.ty) + t_state.tz / (1.0 - rho) + xstr;

        flux_x = 2.0 * psi - flux_x;
        flux_y = 2.0 * psi - flux_y;
        flux_z = (psi - rho * flux_z) / (1.0 - rho);

        return psi;
    }
};
}
}
