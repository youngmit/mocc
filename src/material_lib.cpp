#include "material_lib.hpp"

#include <string>
#include <sstream>
#include <boost/regex.hpp>
#include <iostream>
#include <vector>

#include "material.hpp"
#include "error.hpp"

using std::stringstream;
using std::string;
using std::cout;
using std::cin;

namespace mocc{

MaterialLib::MaterialLib(){
    // do nothing
    return;
}

MaterialLib::MaterialLib(FileScrubber &input){
	string line;
	// Read the first three lines to extract the library header
	// Description
	// Number of groups and number of materials
	{
        // skip the first line
        input.getline();
		stringstream inBuf(input.getline());
		inBuf >> n_grp_;
		if(inBuf.fail()){
			Error("Failed to read number of groups!");
		}
		inBuf >> n_material_;
		if(inBuf.fail()){
			Error("Failed to read number of materials!");
		}
	}

	// Group boundaries
	{
		stringstream inBuf(input.getline());
		for( unsigned int i=0; i<n_grp_; i++ ) {
			double bound;
			inBuf >> bound;
			g_bounds_.push_back(bound);
			if( inBuf.fail() ) {
				Error("Trouble reading group bounds!");
			}
		}
	}
	
	// Read in material data
	for (unsigned int i=0; i<n_material_; i++){
		// Get the name of the material
		line = input.getline();
		
	    boost::regex headExp("^\\s*XSMACRO\\s+([^\\s]+)\\s+([0-9]+)\\s*$");
		boost::smatch results;
		regex_match(line, results, headExp);
        string materialName = results[1].str();
		std::cout << materialName << std::endl;

        // Read in the non-scattering stuff
        VecF abs;
        VecF nuFiss;
        VecF fiss;
        VecF chi;
        for(unsigned int ig=0; ig<n_grp_; ig++){
            stringstream inBuf(input.getline());
            double val;

            inBuf >> val;
            abs.push_back(val);
            inBuf >> val;
            nuFiss.push_back(val);
            inBuf >> val;
            fiss.push_back(val);
            inBuf >> val;
            chi.push_back(val);
            if(inBuf.fail()){
                Error("Trouble reading XS data from library!");
            }
        }

        // Read in the scattering table
        std::vector<VecF> scatTable;
        for(unsigned int ig=0; ig<n_grp_; ig++){
            stringstream inBuf(input.getline());
            VecF scatRow;
            for(unsigned int igg=0; igg<n_grp_; igg++){
                double val;
                inBuf >> val;
                scatRow.push_back(val);
            }
            scatTable.push_back(scatRow);
        }
		
        // produce a Material object and add it to the library
        lib_materials_.insert(std::pair<std::string, Material>(
            materialName.c_str(), 
            Material(abs, nuFiss, fiss, chi, scatTable)));

	}
}

void MaterialLib::assignID(int id, std::string name){
    const Material* mat_p = &lib_materials_.at(name);
    materials_.insert( std::pair<int, const Material*>( id, mat_p) );
    return;
}

};
