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

    private:
    	// Radii of material rings
    	std::vector<mocc::float_t> xs_radii_;
    	// Number of mesh rings
    	int m_nRing;
    	// Radii of actual mesh rings
    	std::vector<mocc::float_t> radii_;
    	// Number of azumuthal subdivisions (for now, for whole pin)
    	std::vector<int> sub_azi_;
    	// Number of radial subdivisions for each material ring
    	std::vector<int> sub_rad_;
    };
}
