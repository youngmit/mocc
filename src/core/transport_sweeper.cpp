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

#include "transport_sweeper.hpp"

#include "pugixml.hpp"

#include <cmath>
#include <iostream>

namespace {
using namespace mocc;
// Locate an ang_quad tag in the XML tree. Start by looking in the current
// node, and consulte parent nodes until an <ang_quad> is found. Return a
// reference to the first <ang_quad> node found. Throw if we get to the
// document root and still dont find one.
const pugi::xml_node find_angquad(const pugi::xml_node &input)
{
    if (input.empty()) {
        throw EXCEPT("Passed node is empty!");
    }

    pugi::xml_node current_node = input;
    while (true) {
        if (!current_node.child("ang_quad").empty()) {
            // We found an <ang_quad>. Return it.
            return current_node.child("ang_quad");
        } else if (!current_node.parent().empty()) {
            // We didnt find an <ang_quad>, but there is a parent node.
            // look there.
            current_node = current_node.parent();
        } else {
            // We reached the end of the line, but still didnt find an
            // <ang_quad>. Fail
            throw EXCEPT(
                "Reached document root without finding an "
                "angular quadrature specification.");
            break;
        }
    }
    return input;
}

SP_XSMesh_t xs_mesh_factory(const CoreMesh &mesh, MeshTreatment treatment)
{
    switch (treatment) {
    case MeshTreatment::TRUE:
    case MeshTreatment::PLANE:
        return std::make_shared<XSMesh>(mesh, treatment);
        break;
    case MeshTreatment::PIN:
    default:
        return std::make_shared<XSMeshHomogenized>(mesh);
    }
}
}

namespace mocc {
TransportSweeper::TransportSweeper(const pugi::xml_node &input,
                                   const CoreMesh &mesh,
                                   MeshTreatment treatment)
    : core_mesh_(&mesh),
      xs_mesh_(xs_mesh_factory(mesh, treatment)),
      n_reg_(mesh.n_reg(treatment)),
      n_group_(xs_mesh_->n_group()),
      groups_(Range(0, n_group_)),
      source_(nullptr),
      flux_(n_reg_, n_group_),
      flux_old_(n_reg_, n_group_),
      vol_(mesh.volumes(treatment)),
      ang_quad_(find_angquad(input)),
      coarse_data_(nullptr),
      n_sweep_(0),
      n_sweep_inner_(0),
      do_incoming_update_(input.attribute("update_incoming").as_bool(true))
{
    return;
}

TransportSweeper::TransportSweeper(const pugi::xml_node &input)
    : source_(nullptr),
      ang_quad_(find_angquad(input)),
      coarse_data_(nullptr),
      n_sweep_(0),
      n_sweep_inner_(0),
      do_incoming_update_(input.attribute("update_incoming").as_bool(true))
{
    return;
}

real_t TransportSweeper::total_fission(bool old) const
{
    real_t tfis      = 0.0;
    const auto &flux = old ? flux_old_ : flux_;
    for (auto &xsr : *xs_mesh_) {
        for (int ig = 0; ig < n_group_; ig++) {
            real_t xsnf = xsr.xsmacnf(ig);
            for (auto &ireg : xsr.reg()) {
                tfis += flux((int)ireg, (int)ig) * vol_[ireg] * xsnf;
            }
        }
    }
    return tfis;
}

void TransportSweeper::calc_fission_source(real_t k,
                                           ArrayB1 &fission_source) const
{
    real_t rkeff   = 1.0 / k;
    fission_source = 0.0;
    for (auto &xsr : *xs_mesh_) {
        const auto &xsnf = xsr.xsmacnf();
        for (int ig = 0; ig < (int)n_group_; ig++) {
            for (auto &ireg : xsr.reg()) {
                fission_source(ireg) += rkeff * xsnf[ig] * flux_(ireg, ig);
            }
        }
    }

    return;
}

/**
 * \todo make this general for all mesh treatments. For now, since this is only
 * used for MoC, we wil just hard-code it to use PLANE treatment.
 */
ArrayB3 TransportSweeper::pin_powers() const
{
    
    assert(n_reg_ == (int)core_mesh_->n_reg(MeshTreatment::PLANE));
    ArrayB3 powers(core_mesh_->subplane().size(), core_mesh_->ny(),
                   core_mesh_->nx());
    powers = 0.0;

    // This isnt the most efficient way to do this, memory-wise, but its
    // quick and simple. Calculate volume x flux x kappa-fission for all
    // flat source regions, then reduce to the pin mesh.
    ArrayB1 fsr_pow(n_reg_);
    fsr_pow  = 0.0;
    int ixsr = 0;
    for (const auto &xsr : *xs_mesh_) {
        for (int ig = 0; ig < n_group_; ig++) {
            for (auto ireg : xsr.reg()) {
                fsr_pow(ireg) += flux_(ireg, ig) * xsr.xsmacf(ig) * vol_[ireg];
            }
        }
        ixsr++;
    }

    int ireg       = 0;
    real_t tot_pow = 0.0;
    int iz         = 0;
    for (int iplane = 0; iplane < (int)core_mesh_->subplane().size();
         ++iplane) {
        auto stt = core_mesh_->begin(iz);
        auto stp = core_mesh_->end(iz);
        int ipin = 0;
        for (auto it = stt; it != stp; ++it) {
            const PinMesh &pm = (*it)->mesh();
            Position pos      = core_mesh_->pin_position(ipin);
            pos.z             = iplane;
            for (int ir = 0; ir < pm.n_reg(); ir++) {
                tot_pow += fsr_pow(ireg);
                powers(iplane, pos.y, pos.x) += fsr_pow(ireg);
                ireg++;
            }
            ipin++;
        }

        iz += core_mesh_->subplane()[iplane];
    }

    // Normalize!
    tot_pow = core_mesh_->n_fuel_2d() / tot_pow;
    for (auto &v : powers) {
        v *= tot_pow;
    }

    return powers;
}

ArrayB2 TransportSweeper::get_pin_flux(MeshTreatment treatment) const
{
    assert(core_mesh_);
    ArrayB2 flux(core_mesh_->n_reg(treatment), n_group_);

    auto flux_it = flux.begin();

    for (int ig = 0; ig < n_group_; ig++) {
        ArrayB1 flux_1g(flux(blitz::Range::all(), ig));
        this->get_pin_flux_1g(ig, flux_1g, treatment);
    }

    return flux;
}

real_t TransportSweeper::flux_residual() const
{
    real_t r    = 0.0;
    auto it     = flux_.begin();
    auto it_old = flux_old_.begin();
    auto end    = flux_.end();
    while (it != end) {
        real_t e = *it - *it_old;
        r += e * e;
        ++it;
        ++it_old;
    }
    return std::sqrt(r);
}
}
