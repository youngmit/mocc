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

#include "coarse_data.hpp"

namespace mocc {
CoarseData::CoarseData(const Mesh &mesh, size_t ngroup)
    : current((int)mesh.n_surf(), ngroup),
      surface_flux((int)mesh.n_surf(), ngroup),
      partial_current((int)mesh.n_surf(), ngroup),
      partial_current_old((int)mesh.n_surf(), ngroup),
      flux((int)mesh.n_pin(), ngroup),
      old_flux((int)mesh.n_pin(), ngroup),
      n_group_(ngroup),
      mesh_(mesh),
      has_data_radial_(false),
      has_data_axial_(false),
      has_old_partial_(false),
      source_("No Data")
{
    current      = 0.0;
    surface_flux = 0.0;
    flux         = 1.0;
    old_flux     = 1.0;
    std::array<real_t, 2> zero = {{0.0, 0.0}};
    partial_current     = zero;
    partial_current_old = zero;

    // assert( current(blitz::Range::all(), 0).isStorageContiguous() );
    // assert( surface_flux(blitz::Range::all(), 0).
    //        isStorageContiguous() );
    // assert( partial_current(blitz::Range::all(), 0).
    //        isStorageContiguous() );

    return;
}

void CoarseData::zero_data(int group, bool zero_partial)
{
    assert(group < n_group_);
    // need this to disambiguate the operator= for partial currents.
    const std::array<real_t, 2> zero = {{0.0, 0.0}};

    current(blitz::Range::all(), group)      = 0.0;
    surface_flux(blitz::Range::all(), group) = 0.0;
    if (zero_partial) {
        partial_current(blitz::Range::all(), group) = zero;
    }
    return;
}

void CoarseData::zero_data_radial(int group, bool zero_partial)
{
    assert(group < n_group_);
    // need this to disambiguate the operator= for partial currents.
    const std::array<real_t, 2> zero = {{0.0, 0.0}};

    auto current_g      = current(blitz::Range::all(), group);
    auto surface_flux_g = surface_flux(blitz::Range::all(), group);
    auto partial_g      = partial_current(blitz::Range::all(), group);

    for (size_t plane = 0; plane < mesh_.nz(); plane++) {
        for (auto surf = mesh_.plane_surf_xy_begin(plane);
             surf != mesh_.plane_surf_end(plane); ++surf) {
            current_g(surf)      = 0.0;
            surface_flux_g(surf) = 0.0;
            if (zero_partial) {
                partial_g(surf) = zero;
            }
        }
    }
    return;
}
}
