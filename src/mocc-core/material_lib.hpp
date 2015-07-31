#pragma once
#include <fstream>
#include <map>
#include <string>
#include "material.hpp"
#include "file_scrubber.hpp"
#include "global_config.hpp"

namespace mocc {
    typedef std::map<unsigned int, const Material*> MaterialMap;


    class MaterialLib{
    public:
    	// Default constructor does nothing
    	MaterialLib();
    	MaterialLib(FileScrubber &input);
    	void assignID(int id, std::string name);

        // Return the number of materials in the library
        unsigned int n_materials() const {
            return n_material_; 
        }

        // Return the map of materials by ID
        const MaterialMap& materials() const {
            return materials_;
        } 

        // Return the number of groups spanned by the library
        unsigned int n_grp() const {
            return n_grp_;
        }

        // Return the group bounds
        const VecF& g_bounds() const {
            return g_bounds_;
        }
    	
    private:
        // Collection of all of the materials specified by the library. Keyed on
        // the string name of the matierial, specified in the first line of the
        // material data blocks
    	std::map<std::string, Material> lib_materials_;
        // Map containing all of the materials that are actually mapped to and
        // ID using <material> tags
        MaterialMap materials_;
        // Number of energy groups for which all materials in the library are
        // defined
    	unsigned int n_grp_;
        // Number of materials that have been mapped to an ID
    	unsigned int n_material_;
        // Number of materials present in the library itself (>= n_material_)
        unsigned int n_material_lib_;
        // Upper bound for each of the energy groups
    	VecF g_bounds_;
        // Descriptive string for the material library
    	std::string m_description;
    };
}
