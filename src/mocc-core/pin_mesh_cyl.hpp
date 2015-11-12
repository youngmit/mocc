#pragma once
#include "global_config.hpp"
#include "pin_mesh_base.hpp"
#include "pugixml.hpp"
#include <vector>

namespace mocc {
    class PinMesh_Cyl : public PinMesh {
    public:
        PinMesh_Cyl(const pugi::xml_node &input);
        ~PinMesh_Cyl();

        int trace( Point2 p1, Point2 p2, int first_reg, VecF &s,
                VecI &reg ) const;

        int find_reg( Point2 p ) const;

        // If i ever get more general with the azimuthal subdivision, i will
        // have to generalize this as well. make sure not to forget.
        size_t n_fsrs( unsigned int xsreg ) const {
            int n = 0;
            if( xsreg<xs_radii_.size() ) {
                n = sub_rad_[xsreg]*sub_azi_[0];
            } else {
                n = sub_azi_[0];
            }
            return n;
        }

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
