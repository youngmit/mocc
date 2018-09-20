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

#include "source.hpp"

#include <algorithm>
#include <iostream>
#include "pugixml.hpp"
#include "util/error.hpp"
#include "core/constants.hpp"

namespace mocc {
Source::Source(int nreg, const XSMesh *xs_mesh, const ArrayB2 &flux)
    : xs_mesh_(xs_mesh),
      n_group_(xs_mesh->n_group()),
      n_reg_(flux.size() / n_group_),
      has_external_(false),
      has_MMS_(false),
      source_1g_(nreg),
      source_1g_with_self_scat_(nreg),
      flux_(flux)
{
    assert(nreg * n_group_ == (int)flux_.size());
    assert(xs_mesh_->n_reg_expanded() == nreg);
    source_1g_.fill(0.0);
    source_1g_with_self_scat_.fill(0.0);
    state_.reset();
    return;
}

void Source::initialize_group(int ig)
{
    if (has_external_) {
        for (int ireg = 0; ireg < (int)n_reg_; ireg++) {
            source_1g_[ireg] = external_source_(ig, ireg);
        }
    }
    else {
        source_1g_.fill(0.0);
    }

    state_.reset();
    return;
}

std::ostream &operator<<(std::ostream &os, const Source &src)
{
    return os;
}

// Multiply the group-independent fission source by \c chi[ig] to get the
// fission source into the current group. If an external source is defines,
// start with that.
void Source::fission(const ArrayB1 &fs, int ig)
{
    assert((int)fs.size() == n_reg_);
    assert(!state_.has_fission);
    assert(!state_.is_scaled);

    for (auto &xsr : *xs_mesh_) {
        real_t xsch = xsr.xsmacch(ig);
        for (const int &ireg : xsr.reg()) {
            source_1g_[ireg] += xsch * fs(ireg);
        }
    }

    state_.has_fission = true;
    return;
}

/**
 * \brief Compute the contribution to the source from inscattering from
 * other groups.
 */
void Source::in_scatter(size_t ig)
{
    assert(!state_.has_inscatter);
    assert(!state_.is_scaled);
    for (auto &xsr : *xs_mesh_) {
        if (xsr.reg().size() == 0) {
            continue;
        }
        const ScatteringRow &scat_row = xsr.xsmacsc().to(ig);
        size_t min_g                  = scat_row.min_g;
        int igg                       = min_g;
        for (auto sc : scat_row) {
            // Dont add a contribution for self-scatter. It might be a
            // good idea to remove self-scatter from the scattering matrix
            // and store it separately. May also benefit the self-scatter
            // routine to have less indirection.
            if (igg != (int)ig) {
                for (auto &ireg : xsr.reg()) {
                    real_t scat_src = sc * flux_((int)ireg, igg);
                    source_1g_[ireg] += scat_src;
                }
            }
            igg++;
        }
    }

    return;
}

void Source::auxiliary(const ArrayB1 &aux)
{
    assert(source_1g_.size() == (int)aux.size());
    assert(!state_.is_scaled);
    for (int i = 0; i < (int)n_reg_; i++) {
        source_1g_[i] += aux(i);
    }
}

void Source::add_external(const pugi::xml_node &input)
{
    if (input.attribute("file").empty()) {
        // Nothing to do here
        return;
    }

    std::string srcfname = input.attribute("file").value();
    H5Node srcfile(srcfname, H5Access::READ);
    srcfile.read("/external_source", external_source_);

    // I converted this from the old HDF5 routines, and the order of these
    // could be wrong. make sure to test or change as necessary
    if (external_source_.extent(0) != (int)n_group_) {
        throw EXCEPT("Wrong group dimensions for source");
    }
    if (external_source_.extent(1) != (int)n_reg_) {
        throw EXCEPT("Wrong regions dimensions for source");
    }

    has_external_ = true;
}

void Source::add_MMS(const pugi::xml_node &input)
{
    if (input.attribute("MMSFile").empty()) {
        // Nothing to do here
        return;
    }

    has_MMS_ = true;

    std::string srcfname = input.attribute("MMSFile").value();

    H5Node MMSh5(srcfname, H5Access::READ);
    // jwg: read in the anisotropic MMS source

    ArrayB3 q_m_j_i;
    MMSh5.read("q_MMS_i_j_m", q_m_j_i);

    auto shape = q_m_j_i.shape();

    int M=shape(0);
    int J=shape(1);
    int I=shape(2);

    ArrayB2 q_ireg_m(I*J,M);

    int ireg=0;
    for (int m=0; m<M; m++){
        ireg=0;
        for (int j=0; j<J; j++){
            for (int i=0; i<I; i++){
                q_ireg_m(ireg,m)=q_m_j_i(m,j,i);
                ireg++;
            }
        }
    }
    MMS_source_.resize(q_ireg_m.shape());
    MMS_source_=q_ireg_m;

    // Read in angular error
    ArrayB2 error_ang_j_i;
    MMSh5.read("error_ang_i_j", error_ang_j_i);

    auto error_ang_shape = error_ang_j_i.shape();
    auto error_ang_size  = error_ang_j_i.size();

    // get the dimensions
    J=error_ang_shape(0);
    I=error_ang_shape(1);

    error_ang_ireg_.resize(error_ang_size,1); // (ireg,ig)
    ireg=0;
    for (int j=0; j<J; j++){
        for (int i=0; i<I; i++){
            error_ang_ireg_(ireg,0)=error_ang_j_i(j,i);
            ireg++;
        }
    }
    // jwg: I left this here because it is strange that the
    // first element assignment generates the following:
    // 2.31779e-310 = 8.88178e-16

    // std::cout << error_ang_j_i(0,0) << std::endl;
    // std::cout << error_ang_ireg_(0,0) << std::endl;
    return;
}
}
