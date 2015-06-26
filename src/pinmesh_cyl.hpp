#pragma once
#include "globalconfig.hpp"
#include "pinmesh_base.hpp"
#include "pugixml.hpp"
#include <vector>

class PinMesh_Cyl : public PinMesh {
public:
	PinMesh_Cyl(const pugi::xml_node &input);
private:
	
	// Number of XS radii
	int m_nRingXS;
	// Radii of material rings
	std::vector<mocc::float_t> m_radiiXS;
	// Number of mesh rings
	int m_nRing;
	// Radii of actual mesh rings
	std::vector<mocc::float_t> m_radii;
	// Number of azumuthal subdivisions (for now, for whole pin)
	std::vector<int> m_subAzi;
	// Number of radial subdivisions for each material ring
	std::vector<int> m_subRad;
	
};