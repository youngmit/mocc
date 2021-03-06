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
#include <iterator>
#include <sstream>
#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"

namespace mocc {
PinMesh_Rect::PinMesh_Rect(const pugi::xml_node &input) : PinMesh(input)
{
    // Parse the number of x and y divisions
    int ndiv_x = input.child("sub_x").text().as_int(0);
    if (ndiv_x < 1) {
        throw EXCEPT("Failed to read valid number of X divisions in rect "
                     "pin mesh.");
    }
    int ndiv_y = input.child("sub_y").text().as_int(0);
    if (ndiv_y < 1) {
        throw EXCEPT("Failed to read valid number of Y divisions in rect "
                     "pin mesh.");
    }

    n_xsreg_ = ndiv_x * ndiv_y;
    n_reg_   = ndiv_x * ndiv_y;

    nx_ = ndiv_x;
    ny_ = ndiv_y;

    real_t dx = pitch_x_ / ndiv_x;
    real_t dy = pitch_y_ / ndiv_y;

    real_t h_pitch_x = 0.5 * pitch_x_;
    real_t h_pitch_y = 0.5 * pitch_y_;

    for (int i = 0; i <= ndiv_x; i++) {
        hx_.push_back(i * dx - h_pitch_x);
    }
    for (int i = 0; i <= ndiv_y; i++) {
        hy_.push_back(i * dy - h_pitch_y);
    }

    // Form lines representing the mesh boundaries
    for (unsigned ix = 1; ix < hx_.size() - 1; ix++) {
        real_t xi = hx_[ix];
        lines_.push_back(Line(Point2(xi, -h_pitch_y), Point2(xi, h_pitch_y)));
    }
    for (unsigned iy = 1; iy < hy_.size() - 1; iy++) {
        real_t yi = hy_[iy];
        lines_.push_back(Line(Point2(-h_pitch_x, yi), Point2(h_pitch_x, yi)));
    }

    // Determine FSR volumes
    areas_ = VecF(n_reg_, dx * dy);

    return;
}

std::pair<real_t, bool> PinMesh_Rect::distance_to_surface(Point2 p,
                                                          Direction dir,
                                                          int &coincident) const
{
    std::pair<real_t, bool> ret;
    int coinc = coincident;

    if ((std::abs(p.x) > 0.5 * pitch_x_) || (std::abs(p.y) > 0.5 * pitch_y_)) {
        ret.first  = 0.0;
        ret.second = true;
        return ret;
    }

    ret.second = false;
    ret.first  = std::numeric_limits<real_t>::max();
    for (const auto &l : lines_) {
        real_t d  = l.distance_to_surface(p, dir, (coincident == l.surf_id));
        if((d < ret.first)) {
            ret.first = d;
            coinc = l.surf_id;
        }
    }
    coincident = coinc;

    return ret;
}

int PinMesh_Rect::trace(Point2 p1, Point2 p2, int first_reg, VecF &s,
                        VecI &reg) const
{
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
int PinMesh_Rect::find_reg(Point2 p) const
{
    // Make sure the point is inside the pin
    if (fabs(p.x) > 0.5 * pitch_x_) {
        return -1;
    }
    if (fabs(p.y) > 0.5 * pitch_y_) {
        return -1;
    }
    
    int ix = std::distance(
        hx_.begin(), std::lower_bound(hx_.cbegin(), hx_.cend(), p.x));
    int iy = std::distance(
        hy_.begin(), std::lower_bound(hy_.cbegin(), hy_.cend(), p.y));
    ix--;
    iy--;
    
    int ireg = nx_ * iy + ix ;

    assert( (ireg >= 0) && (ireg < n_reg_));
    return ireg;
}

int PinMesh_Rect::find_reg(Point2 p, Direction dir) const
{
    // Make sure the point is inside the pin
    if (((p.x < -0.5 * pitch_x_ + REAL_FUZZ) && (dir.ox < 0.0)) ||
        ((p.x > 0.5 * pitch_x_ - REAL_FUZZ) && (dir.ox > 0.0)) ||
        ((p.y < -0.5 * pitch_y_ + REAL_FUZZ) && (dir.oy < 0.0)) ||
        ((p.y > 0.5 * pitch_y_ - REAL_FUZZ) && (dir.oy > 0.0))) {
        return -1;
    }

    int ix = std::distance(
        hx_.begin(), std::lower_bound(hx_.begin(), hx_.end(), p.x, fuzzy_lt));
    if (fp_equiv(p.x, hx_[ix])) {
        ix = (dir.ox > 0.0) ? ix + 1 : ix;
    }
    ix = std::max(0, std::min((int)nx_ - 1, ix-1));
    int iy = std::distance(
        hy_.begin(), std::lower_bound(hy_.begin(), hy_.end(), p.y, fuzzy_lt));
    if (fp_equiv(p.y, hy_[iy])) {
        iy = (dir.oy > 0.0) ? iy + 1 : iy;
    }
    iy = std::max(0, std::min((int)ny_ - 1, iy-1));

    int ireg = nx_ * iy + ix;
    assert(ireg < n_reg_);
    return ireg;
}

void PinMesh_Rect::print(std::ostream &os) const
{
    PinMesh::print(os);
    os << std::endl;
    os << "Type: Rectangular" << std::endl;
    os << "X Divisions: " << std::endl;
    for (const auto &xi : hx_) {
        os << "    " << xi << std::endl;
    }
    os << "Y Divisions: " << std::endl;
    for (const auto &yi : hy_) {
        os << "    " << yi << std::endl;
    }
    return;
}

std::string PinMesh_Rect::draw() const
{
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
