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

#include "core/source.hpp"
#include "core/source_isotropic.hpp"

namespace mocc {
void SourceIsotropic::self_scatter(size_t ig, const ArrayB1 &xstr)
{
    // Take a slice reference for this group's flux
    const ArrayB1 flux_1g = flux_(blitz::Range::all(), ig);
    if (xstr.size() > 0) {
        for (auto &xsr : *xs_mesh_) {
            const ScatteringRow &scat_row = xsr.xsmacsc().to(ig);
            real_t xssc                   = scat_row[ig];
            real_t r_fpi_tr               = 1.0 / (xsr.xsmactr(ig) * FPI);
            for (const int ireg : xsr.reg()) {
                q_[ireg] = (source_1g_[ireg] + flux_1g(ireg) * xssc) * r_fpi_tr;
            }
        }
    } else {
        real_t r_fpi = 1.0 / (FPI);
        for (auto &xsr : *xs_mesh_) {
            const ScatteringRow &scat_row = xsr.xsmacsc().to(ig);
            real_t xssc                   = scat_row[ig];
            for (const int ireg : xsr.reg()) {
                q_[ireg] = (source_1g_[ireg] + flux_1g(ireg) * xssc) * r_fpi;
            }
        }
    }

    // Check to make sure that the source is positive
    bool any = false;
    for (int i = 0; i < q_.size(); i++) {
        if (q_[i] < 0.0) {
            any = true;
        }
    }
    if (any) {
        //  throw EXCEPT("Negative source!");
    }

    return;
}

void SourceIsotropic::self_scatter_for_MMS(size_t ig, const ArrayB1 &xstr)
{
    // Take a slice reference for this group's flux
    const ArrayB1 flux_1g = flux_(blitz::Range::all(), ig);
    if (xstr.size() > 0) {
        for (auto &xsr : *xs_mesh_) {
            const ScatteringRow &scat_row = xsr.xsmacsc().to(ig);
            real_t xssc                   = scat_row[ig];
            real_t r_fpi_tr               = 1.0 / (xsr.xsmactr(ig) * FPI);
            for (const int ireg : xsr.reg()) {
                source_1g_with_self_scat_[ireg]=source_1g_[ireg] + flux_1g(ireg) * xssc;
                // q_1g_[ireg] = source_1g_with_self_scat_[ireg];
                // q_[ireg] = q_1g_[ireg] * r_fpi_tr;
            }
        }
    } else {
        real_t r_fpi = 1.0 / (FPI);
        for (auto &xsr : *xs_mesh_) {
            const ScatteringRow &scat_row = xsr.xsmacsc().to(ig);
            real_t xssc                   = scat_row[ig];
            for (const int ireg : xsr.reg()) {
                source_1g_with_self_scat_[ireg] = source_1g_[ireg] + flux_1g(ireg) * xssc;
                // q_1g_[ireg] = source_1g_with_self_scat_[ireg];
                // q_[ireg] = q_1g_[ireg] * r_fpi;

            }
        }
    }

    // Check to make sure that the source is positive
    bool any = false;
    for (int i = 0; i < q_.size(); i++) {
        if (q_[i] < 0.0) {
            any = true;
        }
    }
    if (any) {
        //  throw EXCEPT("Negative source!");
    }

    return;
}
} // end of namespace mocc
