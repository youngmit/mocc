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

#include "core/error.hpp"
#include "core/geometry/geom.hpp"
#include "core/pin_mesh_base.hpp"
#include "core/pugifwd.hpp"

namespace mocc {
    class PinMesh_Rect: public PinMesh{
    public:
        PinMesh_Rect(const pugi::xml_node &input);

        int trace( Point2 p1, Point2 p2, int first_reg, VecF &s,
                VecI &reg ) const override final;

        int find_reg( Point2 p ) const override final;
        int find_reg( Point2 p, Direction dir ) const override final;

        size_t n_fsrs( unsigned int xsreg ) const {
            return 1;
        }

        virtual std::pair<real_t, Surface> distance_to_surface(Point2 p,
                Direction dir ) const;

        void print( std::ostream &os ) const;

        std::string draw() const;
    private:
        // Vector containing the x divisions, starting at the first internal
        // division, ending at the half-pitch
        VecF hx_;
        // Vector containing the y divisions, defined similarly to the above.
        VecF hy_;
        std::vector<Line> lines_;
    };
}
