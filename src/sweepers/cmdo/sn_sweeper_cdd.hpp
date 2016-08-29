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

#include <memory>
#include "util/files.hpp"
#include "util/force_inline.hpp"
#include "util/pugifwd.hpp"
#include "core/core_mesh.hpp"
#include "core/exponential.hpp"
#include "sn/cell_worker.hpp"
#include "sn/sn_sweeper.hpp"
#include "sn/sn_sweeper_variant.hpp"
#include "correction_data.hpp"

namespace mocc {
namespace cmdo {
using namespace sn;
typedef std::pair<std::unique_ptr<SnSweeper>, std::shared_ptr<CorrectionData>>
    CDDPair_t;

/**
 * \brief Specialization of \ref SnSweeperVariant to use the Corrected Diamond
 * Differenc scheme
 */
template <class Equation>
class SnSweeper_CDD : public SnSweeperVariant<Equation> {
public:
    SnSweeper_CDD(const pugi::xml_node &input, const CoreMesh &mesh)
        : sn::SnSweeperVariant<Equation>(input, mesh), corrections_(nullptr)
    {
        // Look for data to set the angular quadrature
        if (!input.child("data").empty()) {
            std::string fname = input.child("data").attribute("file").value();
            LogScreen << "Reading angular quadrature from file: " << fname
                      << std::endl;
            try {
                H5Node file(fname, H5Access::READ);
                this->ang_quad_ = AngularQuadrature(file);
            } catch (Exception e) {
                throw EXCEPT_E("Failed to read angular quadrature from file",
                               e);
            }

            real_t wsum = 0.0;
            for (auto a : this->ang_quad_) {
                wsum += a.weight;
            }
            std::cout << "Weight sum: " << wsum << std::endl;
        }
        return;
    }

    real_t evaluate_2d(
        real_t &flux_x, real_t &flux_y, real_t q, real_t xstr, int i,
        const typename SnSweeperVariant<Equation>::ThreadState &t_state) const
    {
        int ix    = i % this->mesh_.nx();
        real_t tx = t_state.ox / this->mesh_.dx(ix);

        real_t ax = corrections_->alpha(i, t_state.iang_2d, this->group_,
                                        Normal::X_NORM);
        real_t ay = corrections_->alpha(i, t_state.iang_2d, this->group_,
                                        Normal::Y_NORM);
        real_t b = corrections_->beta(i, t_state.iang_2d, this->group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t psi = q + 2.0 * (tx * flux_x + t_state.ty * flux_y);
        psi /= tx / gx + t_state.ty / gy + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;

        return psi;
    }

    /**
     * \brief Associate the internal reference to correction data.
     *
     * Any existing data will get kicked off. Since this uses
     * std::shared_ptr, if the sweeper has the only reference to any data
     * that gets replaced, we should expect the old data to be destroyed.
     * Usually what we want, but be careful.
     */
    void set_corrections(std::shared_ptr<const CorrectionData> data)
    {
        assert(data);
        corrections_ = data;
    }

protected:
    std::shared_ptr<const CorrectionData> corrections_;
    VecI macroplanes_;
};

/**
 * \brief Specialization of \ref SnSweeper_CDD to use diamond difference in the
 * axial dimension
 */
class SnSweeper_CDD_DD : public SnSweeper_CDD<SnSweeper_CDD_DD> {
public:
    SnSweeper_CDD_DD(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeper_CDD<SnSweeper_CDD_DD>(input, mesh)
    {
        return;
    }

    real_t evaluate_2d(real_t &flux_x, real_t &flux_y, real_t q, real_t xstr,
                       int i, const ThreadState &t_state) const
    {
        return SnSweeper_CDD<SnSweeper_CDD_DD>::evaluate_2d(flux_x, flux_y, q,
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

        real_t psi =
            q + 2.0 * (tx * flux_x + t_state.ty * flux_y + t_state.tz * flux_z);
        psi /= tx / gx + t_state.ty / gy + 2.0 * t_state.tz + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;
        flux_z = 2.0 * psi - flux_z;

        return psi;
    }
};

/**
 * \brief Specialization of \ref SnSweeper_CDD to use forward differencing in
 * the axial dimension
 */
class SnSweeper_CDD_FW : public SnSweeper_CDD<SnSweeper_CDD_FW> {
public:
    SnSweeper_CDD_FW(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeper_CDD<SnSweeper_CDD_FW>(input, mesh)
    {
        return;
    }

    real_t evaluate(real_t &flux_x, real_t &flux_y, real_t &flux_z, real_t q,
                    real_t xstr, int i, const ThreadState &t_state) const
    {
        int ix    = i % mesh_.nx();
        int ia    = t_state.macroplane * plane_size_ + i % plane_size_;
        real_t tx = t_state.ox / mesh_.dx(ix);

        real_t ax =
            corrections_->alpha(ia, t_state.iang_2d, group_, Normal::X_NORM);
        real_t ay =
            corrections_->alpha(ia, t_state.iang_2d, group_, Normal::Y_NORM);
        real_t b = corrections_->beta(ia, t_state.iang_2d, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t psi =
            q + 2.0 * (tx * flux_x + t_state.ty * flux_y) + t_state.tz * flux_z;
        psi /= tx / gx + t_state.ty / gy + t_state.tz + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;
        flux_z = psi;

        return psi;
    }
};

/**
 * \brief Specialization of \ref SnSweeper_CDD to use step characteristics in
 * the axial dimension
 */
class SnSweeper_CDD_SC : public SnSweeper_CDD<SnSweeper_CDD_SC> {
public:
    SnSweeper_CDD_SC(const pugi::xml_node &input, const CoreMesh &mesh)
        : SnSweeper_CDD<SnSweeper_CDD_SC>(input, mesh), exponential_()
    {
        return;
    }

    real_t evaluate(real_t &flux_x, real_t &flux_y, real_t &flux_z, real_t q,
                    real_t xstr, int i, const ThreadState &t_state) const
    {
        int ix    = i % mesh_.nx();
        int ia    = t_state.macroplane * plane_size_ + i % plane_size_;
        real_t tx = t_state.ox / mesh_.dx(ix);

        real_t ax =
            corrections_->alpha(ia, t_state.iang_2d, group_, Normal::X_NORM);
        real_t ay =
            corrections_->alpha(ia, t_state.iang_2d, group_, Normal::Y_NORM);
        real_t b = corrections_->beta(ia, t_state.iang_2d, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t tau    = xstr / t_state.tz;
        real_t rho    = 1.0 / tau - 1.0 / (exponential_.exp(tau) - 1.0);
        real_t rhofac = rho / (1.0 - rho);

        real_t psi = q + 2.0 * (tx * flux_x + t_state.ty * flux_y) +
                     t_state.tz * (rhofac + 1.0) * flux_z;
        psi /= tx / gx + t_state.ty / gy + t_state.tz / (1.0 - rho) + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;
        flux_z = (psi - rho * flux_z) / (1.0 - rho);

        return psi;
    }

private:
    Exponential exponential_;
};
}
}
