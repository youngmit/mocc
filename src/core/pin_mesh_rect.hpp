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

#include "pugixml.hpp"

#include "geom.hpp"
#include "pin_mesh_base.hpp"

namespace mocc {
    class PinMesh_Rect: public PinMesh{
    public:
        PinMesh_Rect(const pugi::xml_node &input);
        int trace( Point2 p1, Point2 p2, int first_reg, VecF &s,
                VecI &reg ) const;

        int find_reg( Point2 p ) const;

        size_t n_fsrs( unsigned int xsreg ) const {
            return 1;
        }

        void print( std::ostream &os ) const;

        std::string draw() const;
    private:
        VecF hx_;
        VecF hy_;
        std::vector<Line> lines_;
    };
}
