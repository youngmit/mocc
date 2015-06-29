#pragma once
#include <fstream>
#include <map>
#include "material.hpp"
#include "file_scrubber.hpp"
#include "global_config.hpp"

namespace mocc {
class MaterialLib{
public:
	// Default constructor does nothing
	MaterialLib();
	MaterialLib(FileScrubber &input);
	void assignID(int id, const char* name);
	
private:
	std::map<const char*, Material> m_materials;
    std::map<int, const char*> m_ids;
	int m_nGrp;
	int m_nMaterial;
	VecF m_gBounds;
	std::string m_description;
};
};
