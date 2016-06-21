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

#include "mesh.hpp"
#include <iostream>
#include <sstream>

namespace mocc {
Mesh::Mesh(size_t n_reg, size_t n_xsreg, VecF &hx, VecF &hy, VecF &hz,
           std::array<Boundary, 6> bc)
    : n_reg_(n_reg),
      n_xsreg_(n_xsreg),
      nx_(hx.size() - 1),
      ny_(hy.size() - 1),
      nz_(hz.size() - 1),
      x_vec_(hx),
      y_vec_(hy),
      z_vec_(hz),
      dx_vec_(hx.size() - 1),
      dy_vec_(hy.size() - 1),
      dz_vec_(hz.size() - 1),
      coarse_vol_(nx_ * ny_ * nz_),
      bounding_box_(Point2(0.0, 0.0), Point2(hx.back(), hy.back())),
      n_surf_plane_((nx_ + 1) * ny_ + (ny_ + 1) * nx_ + nx_ * ny_)
{
    assert(std::is_sorted(x_vec_.begin(), x_vec_.end()));
    assert(std::is_sorted(y_vec_.begin(), y_vec_.end()));
    assert(std::is_sorted(z_vec_.begin(), z_vec_.end()));

    bc_ = bc;

    hx_ = x_vec_.back();
    for (auto xi = x_vec_.cbegin() + 1; xi != x_vec_.cend() - 1; ++xi) {
        lines_.push_back(Line(Point2(*xi, 0.0), Point2(*xi, hx_)));
    }
    hy_ = y_vec_.back();
    for (auto yi = y_vec_.cbegin() + 1; yi != y_vec_.cend() - 1; ++yi) {
        lines_.push_back(Line(Point2(0.0, *yi), Point2(hx_, *yi)));
    }
    hz_ = z_vec_.back();

    bounding_box_ =
        Box(Point2(x_vec_[0], y_vec_[0]), Point2(x_vec_.back(), y_vec_.back()));

    for (size_t i = 1; i < x_vec_.size(); i++) {
        dx_vec_[i - 1] = x_vec_[i] - x_vec_[i - 1];
    }
    for (size_t i = 1; i < y_vec_.size(); i++) {
        dy_vec_[i - 1] = y_vec_[i] - y_vec_[i - 1];
    }
    for (size_t i = 1; i < z_vec_.size(); i++) {
        dz_vec_[i - 1] = z_vec_[i] - z_vec_[i - 1];
    }

    for (size_t i = 0; i < this->n_pin(); i++) {
        auto pos       = this->coarse_position(i);
        coarse_vol_[i] = dx_vec_[pos.x] * dy_vec_[pos.y] * dz_vec_[pos.z];
    }

    assert(nx_ == (int)hx.size() - 1);
    assert(ny_ == (int)hy.size() - 1);
    assert(nz_ == (int)hz.size() - 1);
    this->prepare_surfaces();
    return;
}

void Mesh::prepare_surfaces()
{
    // number of x-normal surfaces in each y row
    int nxsurf = nx_ + 1;
    // number of y-normal surfaces in each x column
    int nysurf = ny_ + 1;
    // number of x- and y-normal surfaces in each plane
    int nxysurf     = (nx_ + 1) * ny_ + (ny_ + 1) * nx_;
    int i           = 0;
    coarse_surf_    = VecI(6 * nx_ * ny_ * nz_, 0);
    int surf_offset = 0;
    for (int iz = 0; iz < nz_; iz++) {
        for (int iy = 0; iy < ny_; iy++) {
            for (int ix = 0; ix < nx_; ix++) {
                Surface surf;
                int cell_offset = i * 6;

                surf = Surface::EAST;
                coarse_surf_[cell_offset + (int)surf] =
                    surf_offset + nxsurf * iy + ix + 1 + nx_ * ny_;

                surf = Surface::WEST;
                coarse_surf_[cell_offset + (int)surf] =
                    surf_offset + nxsurf * iy + ix + nx_ * ny_;

                surf = Surface::NORTH;
                coarse_surf_[cell_offset + (int)surf] =
                    surf_offset + nx_ * ny_ + (nx_ + 1) * ny_ + (ny_ + 1) * ix +
                    iy + 1;

                surf = Surface::SOUTH;
                coarse_surf_[cell_offset + (int)surf] =
                    surf_offset + nxsurf * ny_ + nysurf * ix + iy + nx_ * ny_;

                surf = Surface::BOTTOM;
                coarse_surf_[cell_offset + (int)surf] =
                    surf_offset + nx_ * iy + ix;

                surf = Surface::TOP;
                coarse_surf_[cell_offset + (int)surf] =
                    surf_offset + nx_ * ny_ + nxysurf + nx_ * iy + ix;

                i++;
            }
        }
        surf_offset += nxysurf + nx_ * ny_;
    }
    return;
}

int Mesh::coarse_surf_point(Point2 p, int cell, std::array<int, 2> &s) const
{
    bool on_x = false;
    bool on_y = false;

    int ix = 0;
    for (auto &xi : x_vec_) {
        if (fp_equiv(p.x, xi)) {
            on_x = true;
            break;
        }
        if (xi > p.x) {
            ix--;
            break;
        }
        ix++;
    }

    int iy = 0;
    for (auto &yi : y_vec_) {
        if (fp_equiv(p.y, yi)) {
            on_y = true;
            break;
        }
        if (yi > p.y) {
            iy--;
            break;
        }
        iy++;
    }

    if (on_x && !on_y) {
        // bottom faces (nx_*ny_) PLUS
        // previous rows (nx_+1)*iy PLUS
        // x position (ix)
        // simplified to reduce number of ops
        s[0] = nx_ * (ny_ + iy) + iy + ix;
        s[0] = nx_ * ny_ + (nx_ + 1) * iy + ix;
        return 1;
    }

    if (on_y && !on_x) {
        // bottom faces (nx_*ny_) PLUS
        // all x-normal surfaces (nx_+1)*ny_ PLUS
        // y position (iy)
        // simplified to reduce number of ops, fewer imuls
        s[0] = nx_ * ny_ + nx_ * ny_ + ny_ + iy;
        s[0] = nx_ * ny_ + (nx_ + 1) * ny_ + (ny_ + 1) * ix + iy;
        return 1;
    }

    // If we are on an x-normal and a y-normal face simultaneously, this
    // can be a potential issue for coarse ray data. For each cell in
    // the mesh to have balance, the ray flux leaving/entering a corner
    // must be accounted for, so we can't just say that the flux goes
    // directly into the diagonal neighbor, instead, we must say that it
    // goes into an adjacent neighbor first, then into the diagonal
    // neighbor, even though there are no actual ray segments in the
    // adjacent neighbor. Therefore we may need to return two surface
    // indices. in the event that we are on the border of the geometry, we
    // may only need to return one surface.
    if (on_x && on_y) {
        // For conservation, we need to be consistent with how the ray
        // should cross the corner. It needs to leave the current cell,
        // glance through the adjacent cell, and end up entering on the
        // consistent surface of the diagonal neighbor. We will use the
        // convention that the ray always goes into the x neighbor, grances
        // through, then moves in y to the diagonal neighbor.

        // Find corner of the current cell that we are looking at
        Position cellpos = this->coarse_position(cell);
        Surface corner_x = Surface::INVALID;
        Surface corner_y = Surface::INVALID;

        if (ix == cellpos.x) {
            corner_x = Surface::WEST;
        }
        else if (ix == cellpos.x + 1) {
            corner_x = Surface::EAST;
        }
        if (iy == cellpos.y) {
            corner_y = Surface::SOUTH;
        }
        else if (iy == cellpos.y + 1) {
            corner_y = Surface::NORTH;
        }

        assert(corner_x != Surface::INVALID);
        assert(corner_y != Surface::INVALID);

        Cardinal corner = Cardinal::INVALID;
        if (corner_x == Surface::WEST) {
            if (corner_y == Surface::NORTH) {
                corner = Cardinal::NW;
            }
            else {
                corner = Cardinal::SW;
            }
        }
        else {
            if (corner_y == Surface::NORTH) {
                corner = Cardinal::NE;
            }
            else {
                corner = Cardinal::SE;
            }
        }

        /// So to break down the rules
        /// - on the domain boundary, only return the surface normal to the
        /// boundary. This may need to be re-addressed when we start
        /// thinking about spatial decomposition.
        /// - on the interior, go x normal first, then y normal.

        // Handle the boundary case
        if (ix == 0) {
            int neighbor = this->coarse_neighbor(cell, corner_y);
            switch (corner) {
            case Cardinal::SW:
                s[0] = this->coarse_surf(neighbor, corner_x);
                s[1] = this->coarse_surf(cell, corner_y);
                return 2;
            case Cardinal::NW:
                s[0] = this->coarse_surf(cell, corner_x);
                return 1;
            default:
                assert(false);
            }
        }
        if (ix == nx_) {
            int neighbor = this->coarse_neighbor(cell, corner_y);
            switch (corner) {
            case Cardinal::SE:
                s[0] = this->coarse_surf(neighbor, corner_x);
                s[1] = this->coarse_surf(cell, corner_y);
                return 2;
            case Cardinal::NE:
                s[0] = this->coarse_surf(cell, corner_x);
                return 1;
            default:
                assert(false);
            }
        }
        if (iy == 0) {
            int neighbor = this->coarse_neighbor(cell, corner_x);
            switch (corner) {
            case Cardinal::SW:
                s[0] = this->coarse_surf(neighbor, corner_y);
                s[1] = this->coarse_surf(cell, corner_x);
                return 2;
            case Cardinal::SE:
                s[0] = this->coarse_surf(cell, corner_y);
                return 1;
            default:
                assert(false);
            }
        }
        if (iy == ny_) {
            int neighbor = this->coarse_neighbor(cell, Surface::EAST);
            switch (corner) {
            case Cardinal::NE:
                s[0] = this->coarse_surf(neighbor, Surface::NORTH);
                s[1] = this->coarse_surf(cell, Surface::EAST);
                return 2;
            case Cardinal::NW:
                s[0] = this->coarse_surf(cell, Surface::NORTH);
                return 1;
            default:
                assert(false);
            }
        }

        // So we aren't on a boundary. Handle the interior corner case
        s[0]      = this->coarse_surf(cell, corner_x);
        int neigh = this->coarse_neighbor(cell, corner_x);
        s[1]      = this->coarse_surf(neigh, corner_y);
        return 2;
    }

    return 0;
}

/**
* Given a vector containing two points (which should be on the boundary of
* the mesh), insert points corresponding to intersections of the line formed
* by those points and the interfaces of all of the cells in the Mesh. The
* points are added to the passed vector and sorted.
*/
void Mesh::trace(std::vector<Point2> &ps) const
{
    assert(ps.size() == 2);

    Point2 p1 = ps[0];
    Point2 p2 = ps[1];
    assert(p2.y > p1.y);

    Line l(p1, p2);

    for (auto &li : lines_) {
        Point2 intersection;
        if (Intersect(li, l, intersection) == 1) {
            ps.push_back(intersection);
        }
    }

    // Since the original points passed in are already on the domain boundary,
    // skip intersections with the bounding box.

    // Sort the points and remove duplicates
    std::sort(ps.begin(), ps.end());
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());

    return;
}

/**
 * \brief Return a \ref Position corresponding to a \ref Point2.
 *
 * The returned \ref Position is always in the zero-th plane.
 * It is possible to produce an invalid \ref Position, if the supplied point
 * does not lie in the domain.
 */
int Mesh::coarse_cell_point(Point2 p) const
{
    auto ix = std::lower_bound(x_vec_.begin(), x_vec_.end(), p.x) -
              x_vec_.begin() - 1;
    auto iy =
        std::distance(y_vec_.begin(),
                      std::lower_bound(y_vec_.begin(), y_vec_.end(), p.y)) -
        1;

    return this->coarse_cell(Position(ix, iy, 0));
}

/**
 * \brief Return a \ref Position corresponding to a \ref Point3.
 *
 * It is possible to produce an invalid \ref Position, if the supplied point
 * does not lie in the domain.
 */
int Mesh::coarse_cell_point(Point3 p) const
{
    auto ix = std::lower_bound(x_vec_.begin(), x_vec_.end(), p.x) -
              x_vec_.begin() - 1;
    auto iy =
        std::distance(y_vec_.begin(),
                      std::lower_bound(y_vec_.begin(), y_vec_.end(), p.y)) -
        1;
    auto iz =
        std::distance(z_vec_.begin(),
                      std::lower_bound(z_vec_.begin(), z_vec_.end(), p.z)) -
        1;

    /// \todo add z support
    return this->coarse_cell(Position(ix, iy, iz));
}

size_t Mesh::coarse_norm_point(Point2 p, int octant, Surface (&s)[2]) const
{
    assert(octant < 5);

    bool on_x = false;
    bool on_y = false;

    int ix = 0;
    for (auto &xi : x_vec_) {
        if (fp_equiv(p.x, xi)) {
            on_x = true;
            break;
        }
        if (xi > p.x) {
            ix--;
            break;
        }
        ix++;
    }

    int iy = 0;
    for (auto &yi : y_vec_) {
        if (fp_equiv(p.y, yi)) {
            on_y = true;
            break;
        }
        if (yi > p.y) {
            iy--;
            break;
        }
        iy++;
    }

    // Return super early if we arent even on a face
    if (!on_x && !on_y) {
        return 0;
    }

    // Return early if we have a clean intersection
    if (on_x != on_y) {
        // Upwind domain boundaries are a little different
        if (on_x) {
            if ((ix == 0) && ((octant == 1) || (octant == 4))) {
                s[0] = Surface::WEST;
                return 1;
            }
            if ((ix == nx_) && ((octant == 2) || (octant == 3))) {
                s[0] = Surface::EAST;
                return 1;
            }
        }
        if (on_y) {
            if ((iy == 0) && ((octant == 1) || (octant == 2))) {
                s[0] = Surface::SOUTH;
                return 1;
            }
            if ((iy == ny_) && ((octant == 3) || (octant == 4))) {
                s[0] = Surface::NORTH;
                return 1;
            }
        }

        // Interior crossings (yay, some sanity!)
        if (on_x) {
            if ((octant == 1) || (octant == 4)) {
                s[0] = Surface::EAST;
            }
            if ((octant == 2) || (octant == 3)) {
                s[0] = Surface::WEST;
            }
        }
        else {
            if ((octant == 1) || (octant == 2)) {
                s[0] = Surface::NORTH;
            }
            if ((octant == 3) || (octant == 4)) {
                s[0] = Surface::SOUTH;
            }
        }
        return 1;
    }

    // Interior corners are a little different than boundary corners
    if ((ix > 0) && (ix < nx_) && (iy > 0) && (iy < ny_)) {
        switch (octant) {
        case 1:
            s[0] = Surface::EAST;
            s[1] = Surface::NORTH;
            break;
        case 2:
            s[0] = Surface::WEST;
            s[1] = Surface::NORTH;
            break;
        case 3:
            s[0] = Surface::WEST;
            s[1] = Surface::SOUTH;
            break;
        case 4:
            s[0] = Surface::EAST;
            s[1] = Surface::SOUTH;
            break;
        }
        return 2;
    }

    // If we make it this far, our point is sitting on a boundary corner
    if (ix == 0) {
        switch (octant) {
        case 1:
            s[0] = Surface::WEST;
            return 1;
        case 2:
            s[0] = Surface::NORTH;
            s[1] = Surface::WEST;
            return 2;
        case 3:
            s[0] = Surface::SOUTH;
            s[1] = Surface::WEST;
            return 2;
        case 4:
            s[0] = Surface::WEST;
            return 1;
        }
    }

    if (ix == nx_) {
        switch (octant) {
        case 1:
            s[0] = Surface::NORTH;
            s[1] = Surface::EAST;
            return 2;
        case 2:
            s[0] = Surface::EAST;
            return 1;
        case 3:
            s[0] = Surface::EAST;
            return 1;
        case 4:
            s[0] = Surface::SOUTH;
            s[1] = Surface::EAST;
            return 2;
        }
    }

    if (iy == 0) {
        switch (octant) {
        case 1:
            s[0] = Surface::SOUTH;
            return 1;
        case 2:
            s[0] = Surface::SOUTH;
            return 1;
        case 3:
            s[0] = Surface::WEST;
            s[1] = Surface::SOUTH;
            return 2;
        case 4:
            s[0] = Surface::EAST;
            s[1] = Surface::SOUTH;
            return 2;
        }
    }

    if (iy == ny_) {
        switch (octant) {
        case 1:
            s[0] = Surface::EAST;
            s[1] = Surface::NORTH;
            return 2;
        case 2:
            s[0] = Surface::WEST;
            s[1] = Surface::NORTH;
            return 2;
        case 3:
            s[0] = Surface::NORTH;
            return 1;
        case 4:
            s[0] = Surface::NORTH;
            return 1;
        }
    }

    return 0;
}

int Mesh::coarse_boundary_cell(Point2 p, int octant) const
{
    assert(octant < 5);

    bool on_x = false;
    bool on_y = false;

    int ix = 0;
    for (auto &xi : x_vec_) {
        if (fp_equiv(p.x, xi)) {
            on_x = true;
            break;
        }
        if (xi > p.x) {
            ix--;
            break;
        }
        ix++;
    }

    int iy = 0;
    for (auto &yi : y_vec_) {
        if (fp_equiv(p.y, yi)) {
            on_y = true;
            break;
        }
        if (yi > p.y) {
            iy--;
            break;
        }
        iy++;
    }

    assert(ix <= nx_);
    assert(iy <= ny_);

    // Determine which boundary we are on
    if (fp_equiv_abs(p.x, 0.0)) {
        // On the west side
        assert((octant == 1) || (octant == 4));
        if (octant == 1) {
            // Good to leave as-is
        }
        else {
            if (on_y) {
                // As per convention, bump down one in y
                iy--;
            }
        }
    }
    else if (fp_equiv(p.x, hx_)) {
        // On the east side
        ix--;
        assert((octant == 2) || (octant == 3));
        if (octant == 2) {
            // Good to leave as-is
        }
        else {
            if (on_y) {
                // As per convention, bump down one in y
                iy--;
            }
        }
    }
    else if (fp_equiv_abs(p.y, 0.0)) {
        // On the south side
        assert((octant == 1) || (octant == 2));
        if (octant == 1) {
            // Good to leave as-is
        }
        else {
            if (on_x) {
                // As per convention, bump down one in x
                ix--;
            }
        }
    }
    else if (fp_equiv(p.y, hy_)) {
        // On the north side
        iy--;
        assert((octant == 3) || (octant == 4));
        if (octant == 3) {
            if (on_x) {
                // As per convention, bump down one in x
                ix--;
            }
        }
        else {
            // Good to leave as-is
        }
    }
    else {
        assert(false);
    }

    size_t cell = this->coarse_cell(Position(ix, iy, 0));

    if ((cell < 0) || (cell >= this->n_pin())) {
        std::stringstream msg;
        msg << cell << " " << ix << " " << iy << std::endl;
        msg << p.x << " " << p.y;
        throw EXCEPT(msg.str());
    }

    return cell;
}
}
