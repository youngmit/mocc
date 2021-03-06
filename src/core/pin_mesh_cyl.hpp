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

#include <string>
#include <vector>
#include "util/global_config.hpp"
#include "util/pugifwd.hpp"
#include "core/pin_mesh_base.hpp"

namespace mocc {
class PinMesh_Cyl : public PinMesh {
public:
    PinMesh_Cyl(const pugi::xml_node &input);
    ~PinMesh_Cyl();

    int trace(Point2 p1, Point2 p2, int first_reg, VecF &s,
              VecI &reg) const final override;

    int find_reg(Point2 p) const final override;
    int find_reg(Point2 p, Direction dir) const final override;

    // If i ever get more general with the azimuthal subdivision, i will
    // have to generalize this as well.
    size_t n_fsrs(unsigned int xsreg) const override
    {
        int n = 0;
        if (xsreg < xs_radii_.size()) {
            n = sub_rad_[xsreg] * sub_azi_[0];
        }
        else {
            n = sub_azi_[0];
        }
        return n;
    }

    std::pair<real_t, bool> distance_to_surface(Point2 p, Direction dir,
                                                int &coincident) const override;

    void print(std::ostream &os) const override;

    std::string draw() const override;

private:
    // Radii of material rings
    std::vector<real_t> xs_radii_;
    // Radii of actual mesh rings
    std::vector<real_t> radii_;
    // Vector of circle objects.
    std::vector<Circle> circles_;
    // Vector of line objects
    std::vector<Line> lines_;
    // Number of azumuthal subdivisions (for now, for whole pin)
    VecI sub_azi_;
    // Number of radial subdivisions for each material ring
    VecI sub_rad_;
};
}
