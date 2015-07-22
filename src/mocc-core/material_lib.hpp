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
    	
    private:
    	std::map<std::string, Material> lib_materials_;
        MaterialMap materials_;
    	unsigned int n_grp_;
    	unsigned int n_material_;
    	VecF g_bounds_;
    	std::string m_description;
    };
}
