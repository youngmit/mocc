#pragma once
#include "pin_mesh_base.hpp"

#include <vector>

#include "pugixml.hpp"

#include "geom.hpp"



namespace mocc {
    class PinMesh_Rect: public PinMesh{
    public:
    	PinMesh_Rect(const pugi::xml_node &input);
        int trace( Point2 p1, Point2 p2, int first_reg, VecF &s, 
                VecI &reg ) const;

        int find_reg( Point2 p ) const;

        int n_fsrs( unsigned int xsreg ) const {
            return 1;
        }
    private:
        VecF hx_;
        VecF hy_;
        std::vector<Line> lines_;
    };
}
