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

#include "pin_mesh_rect.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include "pugixml.hpp"

#include "error.hpp"

namespace mocc {
PinMesh_Rect::PinMesh_Rect(const pugi::xml_node &input) : PinMesh(input) {
    // Parse the number of x and y divisions
    int ndiv_x = input.child("sub_x").text().as_int(0);
    if (ndiv_x < 1) {
        throw EXCEPT(
            "Failed to read valid number of X divisions in rect "
            "pin mesh.");
    }
    int ndiv_y = input.child("sub_y").text().as_int(0);
    if (ndiv_y < 1) {
        throw EXCEPT(
            "Failed to read valid number of Y divisions in rect "
            "pin mesh.");
    }

    n_xsreg_ = ndiv_x * ndiv_y;
    n_reg_ = ndiv_x * ndiv_y;

    real_t dx = pitch_x_ / ndiv_x;
    real_t dy = pitch_y_ / ndiv_y;

    real_t h_pitch_x = 0.5 * pitch_x_;
    real_t h_pitch_y = 0.5 * pitch_y_;

    for (int i = 1; i < ndiv_x; i++) {
        hx_.push_back(i * dx - h_pitch_x);
    }
    for (int i = 1; i < ndiv_y; i++) {
        hy_.push_back(i * dy - h_pitch_y);
    }

    // Form lines representing the mesh boundaries
    for (auto &xi : hx_) {
        lines_.push_back(Line(Point2(xi, -h_pitch_y), Point2(xi, h_pitch_y)));
    }
    for (auto &yi : hy_) {
        lines_.push_back(Line(Point2(-h_pitch_x, yi), Point2(h_pitch_x, yi)));
    }

    // Determine FSR volumes
    vol_ = VecF(n_reg_, dx * dy);

    return;
}

std::pair<real_t, Surface> PinMesh_Rect::distance_to_surface(Point2 p,
                                                          Direction dir) const {
    std::pair<real_t, Surface> ret;
    ret.second = Surface::INTERNAL;
    ret.first = std::numeric_limits<real_t>::max();
    for (const auto &l : lines_) {
        real_t d = l.distance_to_surface(p, dir);
        ret.first = (d < ret.first) ? d : ret.first;
    }

    Box bound(Point2(-0.5 * pitch_x_, -0.5 * pitch_y_),
              Point2(0.5 * pitch_x_, 0.5 * pitch_y_));
    auto bound_d = bound.distance_to_surface(p, dir);
    if (bound_d.first < ret.first) {
        ret.first = bound_d.first;
        ret.second = bound_d.second;
    }

    return ret;
}

int PinMesh_Rect::trace(Point2 p1, Point2 p2, int first_reg, VecF &s,
                        VecI &reg) const {
    // Make a vector to store the collision points and add the start and end
    // points to it
    std::vector<Point2> ps;
    ps.push_back(p1);
    ps.push_back(p2);

    // Create a line object for the input points to test for collisions with
    Line l(p1, p2);

    // Test for collisions with all of the lines in the mesh
    for (auto &li : lines_) {
        Point2 p;
        int ret = Intersect(li, l, p);
        if (ret == 1) {
            ps.push_back(p);
        }
    }

    // Sort the points
    std::sort(ps.begin(), ps.end());
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());

    // Determine segment lengths and region indices
    for (unsigned int ip = 1; ip < ps.size(); ip++) {
        s.push_back(ps[ip].distance(ps[ip - 1]));
        reg.push_back(this->find_reg(Midpoint(ps[ip], ps[ip - 1])) + first_reg);
    }

    return ps.size() - 1;
}

/**
 * In the rectangular mesh, the indices are ordered naturally. The first
 * region is in the lower left, the last in the upper right, proceeding
 * first in the x direction, then in the y. nothing too fancy.
 */
int PinMesh_Rect::find_reg(Point2 p) const {
    // Make sure the point is inside the pin
    if (fabs(p.x) > 0.5 * pitch_x_) {
        return -1;
    }
    if (fabs(p.y) > 0.5 * pitch_y_) {
        return -1;
    }

    int ix;
    for (ix = 0; ix < (int)hx_.size(); ix++) {
        if (hx_[ix] > p.x) {
            break;
        }
    }
    int iy;
    for (iy = 0; iy < (int)hy_.size(); iy++) {
        if (hy_[iy] > p.y) {
            break;
        }
    }

    int ireg = (hx_.size() + 1) * iy + ix;
    assert(ireg < n_reg_);
    return ireg;
}

void PinMesh_Rect::print(std::ostream &os) const {
    PinMesh::print(os);
    os << std::endl;
    os << "Type: Rectangular" << std::endl;
    os << "X Divisions: " << std::endl;
    for( const auto &xi: hx_ ) {
        os << "    " << xi << std::endl;
    }
    os << "Y Divisions: " << std::endl;
    for( const auto &yi: hy_ ) {
        os << "    " << yi << std::endl;
    }
    return;
}

int PinMesh_Rect::find_reg(Point2 p, Direction dir) const {
    // Make sure the point is inside the pin
    if (fabs(p.x) > 0.5 * pitch_x_) {
        return -1;
    }
    if (fabs(p.y) > 0.5 * pitch_y_) {
        return -1;
    }

    int ix = std::distance(
        hx_.begin(), std::lower_bound(hx_.begin(), hx_.end(), p.x, fuzzy_lt));
    if (fp_equiv_abs(p.x, hx_[ix])) {
        ix = (dir.ox > 0.0) ? ix + 1 : ix;
    }
    int iy = std::distance(
        hy_.begin(), std::lower_bound(hy_.begin(), hy_.end(), p.y, fuzzy_lt));
    if (fp_equiv_abs(p.y, hy_[iy])) {
        iy = (dir.oy > 0.0) ? iy + 1 : iy;
    }
    int ireg = (hx_.size() + 1) * iy + ix;
    assert(ireg < n_reg_);
    return ireg;
}

std::string PinMesh_Rect::draw() const {
    std::stringstream buf;

    for (auto l : lines_) {
        buf << "ctx.move_to(" << l.p1.x << ", " << l.p1.y << ")" << std::endl;
        buf << "ctx.line_to(" << l.p2.x << ", " << l.p2.y << ")" << std::endl;
        buf << "ctx.close_path()" << std::endl;
    }
    buf << "ctx.stroke()";

    return buf.str();
}
}
