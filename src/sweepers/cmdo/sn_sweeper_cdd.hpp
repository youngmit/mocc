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
 * An extension of \ref sn::CellWorker to propagate flux through an
 * orthogonal mesh region with the corrected diamond difference (CDD)
 * scheme. This class is still virtual, as the
 * sn::CellWorker::evaluate() method can be tailored for different axial
 * treatments.
 */
class CellWorker_CDD : public sn::CellWorker {
public:
    CellWorker_CDD(const Mesh &mesh, const AngularQuadrature &ang_quad)
        : CellWorker(mesh, ang_quad), ang_quad_(ang_quad), corrections_(nullptr)
    {
        return;
    }

    inline virtual void set_z(int iz) override final
    {
        CellWorker::set_z(iz);
        correction_z_ = mesh_.macroplane_index()[iz];
    }

    void set_group(int group) override final
    {
        group_ = group;
    }

    MOCC_FORCE_INLINE void set_angle(int iang, Angle angle) override final
    {
        sn::CellWorker::set_angle(iang, angle);
        iang_alpha_ = iang % (ang_quad_.ndir() / 2);
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
        corrections_ = data;
    }

    /**
     * \copydoc sn::CellWorker::evaluate_2d()
     *
     * Since the variants of the CDD worker are all for different axial
     * treatments, the 2-D version of \ref sn::CellWorker::evaluate() can
     * live here.
     */
    MOCC_FORCE_INLINE real_t evaluate_2d(real_t &flux_x, real_t &flux_y,
                                         real_t q, real_t xstr, int i) override
    {
        int ix    = i % mesh_.nx();
        real_t tx = ox_ / mesh_.dx(ix);

        real_t ax = corrections_->alpha(i, iang_alpha_, group_, Normal::X_NORM);
        real_t ay = corrections_->alpha(i, iang_alpha_, group_, Normal::Y_NORM);
        real_t b  = corrections_->beta(i, iang_alpha_, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t psi = q + 2.0 * (tx * flux_x + ty_ * flux_y);
        psi /= tx / gx + ty_ / gy + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;

        return psi;
    }

protected:
    const AngularQuadrature &ang_quad_;

    std::shared_ptr<const CorrectionData> corrections_;

    int iang_alpha_;
    int correction_z_;

    int group_;
};

/**
 * An extension of \ref CellWorker_CDD to propagate flux through an
 * orthogonal mesh region with the corrected diamond difference (CDD)
 * in X and Y, with diamond difference in Z.
 */
class CellWorker_CDD_DD : public CellWorker_CDD {
public:
    CellWorker_CDD_DD(const Mesh &mesh, const AngularQuadrature &ang_quad)
        : CellWorker_CDD(mesh, ang_quad)
    {
        return;
    }

    MOCC_FORCE_INLINE real_t evaluate(real_t &flux_x, real_t &flux_y,
                                      real_t &flux_z, real_t q, real_t xstr,
                                      int i) override final
    {
        int ix    = i % mesh_.nx();
        int ia    = correction_z_ * plane_size_ + i % plane_size_;
        real_t tx = ox_ / mesh_.dx(ix);

        real_t ax =
            corrections_->alpha(ia, iang_alpha_, group_, Normal::X_NORM);
        real_t ay =
            corrections_->alpha(ia, iang_alpha_, group_, Normal::Y_NORM);
        real_t b = corrections_->beta(ia, iang_alpha_, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t psi = q + 2.0 * (tx * flux_x + ty_ * flux_y + tz_ * flux_z);
        psi /= tx / gx + ty_ / gy + 2.0 * tz_ + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;
        flux_z = 2.0 * psi - flux_z;

        return psi;
    }
};

/**
 * A variant of \ref CellWorker_CDD to propagate flux through an orthogonal
 * mesh region with the corrected diamond difference (CDD) scheme in X and
 * Y, with FW difference in Z.
 */
class CellWorker_CDD_FW : public CellWorker_CDD {
public:
    CellWorker_CDD_FW(const Mesh &mesh, const AngularQuadrature &ang_quad)
        : CellWorker_CDD(mesh, ang_quad)
    {
        return;
    }

    MOCC_FORCE_INLINE real_t evaluate(real_t &flux_x, real_t &flux_y,
                                      real_t &flux_z, real_t q, real_t xstr,
                                      int i) override final
    {
        int ix    = i % mesh_.nx();
        int ia    = correction_z_ * plane_size_ + i % plane_size_;
        real_t tx = ox_ / mesh_.dx(ix);

        real_t ax =
            corrections_->alpha(ia, iang_alpha_, group_, Normal::X_NORM);
        real_t ay =
            corrections_->alpha(ia, iang_alpha_, group_, Normal::Y_NORM);
        real_t b = corrections_->beta(ia, iang_alpha_, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t psi = q + 2.0 * (tx * flux_x + ty_ * flux_y) + tz_ * flux_z;
        psi /= tx / gx + ty_ / gy + tz_ + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;
        flux_z = psi;

        return psi;
    }
};

/**
 * A variant of \ref CellWorker_CDD to propagate flux through an orthogonal
 * mesh region with the corrected diamond difference (CDD) scheme in X and
 * Y, with step characteristics in Z.
 */
class CellWorker_CDD_SC : public CellWorker_CDD {
public:
    CellWorker_CDD_SC(const Mesh &mesh, const AngularQuadrature &ang_quad)
        : CellWorker_CDD(mesh, ang_quad), exponential_()
    {
        return;
    }

    MOCC_FORCE_INLINE real_t evaluate(real_t &flux_x, real_t &flux_y,
                                      real_t &flux_z, real_t q, real_t xstr,
                                      int i) override final
    {
        int ix    = i % mesh_.nx();
        int ia    = correction_z_ * plane_size_ + i % plane_size_;
        real_t tx = ox_ / mesh_.dx(ix);

        real_t ax =
            corrections_->alpha(ia, iang_alpha_, group_, Normal::X_NORM);
        real_t ay =
            corrections_->alpha(ia, iang_alpha_, group_, Normal::Y_NORM);
        real_t b = corrections_->beta(ia, iang_alpha_, group_);

        real_t gx = ax * b;
        real_t gy = ay * b;

        real_t tau    = xstr / tz_;
        real_t rho    = 1.0 / tau - 1.0 / (exponential_.exp(tau) - 1.0);
        real_t rhofac = rho / (1.0 - rho);

        real_t psi = q + 2.0 * (tx * flux_x + ty_ * flux_y) +
                     tz_ * (rhofac + 1.0) * flux_z;
        psi /= tx / gx + ty_ / gy + tz_ / (1.0 - rho) + xstr;

        flux_x = (psi - gx * flux_x) / gx;
        flux_y = (psi - gy * flux_y) / gy;
        flux_z = (psi - rho * flux_z) / (1.0 - rho);

        return psi;
    }

private:
    Exponential exponential_;
};

template <class T> class SnSweeper_CDD : public sn::SnSweeperVariant<T> {
public:
    SnSweeper_CDD(const pugi::xml_node &input, const CoreMesh &mesh)
        : sn::SnSweeperVariant<T>(input, mesh), correction_data_(nullptr)
    {
        // Look for data to set the angular quadrature
        if (!input.child("data").empty()) {
            std::string fname = input.child("data").attribute("file").value();
            H5Node file(fname, H5Access::READ);
            this->ang_quad_ = AngularQuadrature(file);
            real_t wsum     = 0.0;
            for (auto a : this->ang_quad_) {
                wsum += a.weight;
            }
            std::cout << "Weight sum: " << wsum << std::endl;
        }
        return;
    }

    void set_corrections(std::shared_ptr<const CorrectionData> data)
    {
        correction_data_ = data;
        this->cell_worker_.set_corrections(data);
    }

private:
    std::shared_ptr<const CorrectionData> correction_data_;
};
}
}
