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

#include "correction_worker.hpp"

#include <cmath>
#include <iomanip>

using std::cout;
using std::endl;

namespace mocc {
namespace cmdo {
/**
 * Quick note on indexing: All of the buffer arrays storing the flux and
 * cross-section sums are sized to a single plane's worth of cells and
 * surfaces, since the angle loop is inside the plane loop, and this
 * routine is called for each angle, for each plane. Therefore, we use
 * an offset to store data into the mesh-global correction factor
 * storage, but none to access the single-plane buffers.
 */
void CurrentCorrections::calculate_corrections(size_t ang, size_t group)
{
    const int FW = 0;
    const int BW = 1;

    const int XL = 0;
    const int XR = 1;
    const int YL = 2;
    const int YR = 3;
    int iang1    = ang;
    int iang2    = ang_quad_.reverse(ang);
    real_t ox    = ang_quad_[ang].ox;

    Surface surfs[2][4];
    // We know that all of our moc angles are positive in the y
    // direction
    surfs[FW][YL] = Surface::SOUTH;
    surfs[FW][YR] = Surface::NORTH;
    surfs[BW][YL] = Surface::NORTH;
    surfs[BW][YR] = Surface::SOUTH;

    if (ox > 0.0) {
        surfs[FW][XL] = Surface::WEST;
        surfs[FW][XR] = Surface::EAST;
        surfs[BW][XL] = Surface::EAST;
        surfs[BW][XR] = Surface::WEST;
    }
    else {
        surfs[FW][XL] = Surface::EAST;
        surfs[FW][XR] = Surface::WEST;
        surfs[BW][XL] = Surface::WEST;
        surfs[BW][XR] = Surface::EAST;
    }

    // See the Surface Normalization page
    // Note that the sin and cos are flipped from what we have in \ref
    // RayData::RayData(). This is because "x spacing" applies to the y-normal
    // faces and vice versa.  we are currently using a ray count-based
    // normalization. so the above doesnt really apply
    real_t area[2] = {std::abs(rays_.spacing(ang) / cos(ang_quad_[ang].alpha)),
                      std::abs(rays_.spacing(ang) / sin(ang_quad_[ang].alpha))};
    for (size_t ic = 0; ic < mesh_->n_cell_plane(); ic++) {
        int icc = ic + cell_offset_;

        auto pos = mesh_->coarse_position(ic);

        real_t area_x = area[0] / mesh_->pin_dx()[pos.x];
        real_t area_y = area[1] / mesh_->pin_dy()[pos.y];

        real_t xstr = sn_xs_mesh_[icc].xsmactr()[group];

        // FW direction
        {
            real_t psi_xl =
                surf_sum_(mesh_->coarse_surf(ic, surfs[FW][XL]) * 2 + 0) *
                area_x;
            real_t psi_xr =
                surf_sum_(mesh_->coarse_surf(ic, surfs[FW][XR]) * 2 + 0) *
                area_x;
            real_t psi_yl =
                surf_sum_(mesh_->coarse_surf(ic, surfs[FW][YL]) * 2 + 0) *
                area_y;
            real_t psi_yr =
                surf_sum_(mesh_->coarse_surf(ic, surfs[FW][YR]) * 2 + 0) *
                area_y;
            real_t ax = vol_sum_(ic * 2 + 0) / (psi_xl + psi_xr);
            real_t ay = vol_sum_(ic * 2 + 0) / (psi_yl + psi_yr);

            real_t b = sigt_sum_(ic * 2 + 0) / xstr;

            real_t e =
                ax - corrections_->alpha(icc, iang1, group, Normal::X_NORM);
            residual_[0] += e * e;
            e = ay - corrections_->alpha(icc, iang1, group, Normal::Y_NORM);
            residual_[1] += e * e;
            e = b - corrections_->beta(icc, iang1, group);
            residual_[2] += e * e;

            corrections_->alpha(icc, iang1, group, Normal::X_NORM) = ax;
            corrections_->alpha(icc, iang1, group, Normal::Y_NORM) = ay;
            corrections_->beta(icc, iang1, group) = b;
        }

        // BW direction
        {
            real_t psi_xl =
                surf_sum_(mesh_->coarse_surf(ic, surfs[BW][XL]) * 2 + 1) *
                area_x;
            real_t psi_xr =
                surf_sum_(mesh_->coarse_surf(ic, surfs[BW][XR]) * 2 + 1) *
                area_x;
            real_t psi_yl =
                surf_sum_(mesh_->coarse_surf(ic, surfs[BW][YL]) * 2 + 1) *
                area_y;
            real_t psi_yr =
                surf_sum_(mesh_->coarse_surf(ic, surfs[BW][YR]) * 2 + 1) *
                area_y;

            real_t ax = vol_sum_(ic * 2 + 1) / (psi_xl + psi_xr);
            real_t ay = vol_sum_(ic * 2 + 1) / (psi_yl + psi_yr);

            real_t b = sigt_sum_(ic * 2 + 1) / xstr;

            real_t e =
                ax - corrections_->alpha(icc, iang2, group, Normal::X_NORM);
            residual_[0] += e * e;
            e = ay - corrections_->alpha(icc, iang2, group, Normal::Y_NORM);
            residual_[1] += e * e;
            e = b - corrections_->beta(icc, iang2, group);
            residual_[2] += e * e;

            corrections_->alpha(icc, iang2, group, Normal::X_NORM) = ax;
            corrections_->alpha(icc, iang2, group, Normal::Y_NORM) = ay;

            corrections_->beta(icc, iang2, group) = b;
        }
    }

    return;
}
}
}
