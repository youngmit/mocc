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
		inBuf >> m_nGrp;
		if(inBuf.fail()){
			Error("Failed to read number of groups!");
		}
		inBuf >> m_nMaterial;
		if(inBuf.fail()){
			Error("Failed to read number of materials!");
		}
	}

	// Group boundaries
	{
		stringstream inBuf(input.getline());
		for (int i=0; i<m_nGrp; i++){
			double bound;
			inBuf >> bound;
			m_gBounds.push_back(bound);
			if(inBuf.fail()){
				Error("Trouble reading group bounds!");
			}
		}
	}
	
	// Read in material data
	for (int i=0; i<m_nMaterial; i++){
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
        for(int ig=0; ig<m_nGrp; ig++){
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
        for(int ig=0; ig<m_nGrp; ig++){
            stringstream inBuf(input.getline());
            VecF scatRow;
            for(int igg=0; igg<m_nGrp; igg++){
                double val;
                inBuf >> val;
                scatRow.push_back(val);
            }
            scatTable.push_back(scatRow);
        }
		
        // produce a Material object and add it to the library
        m_materials.insert(std::pair<const char*, Material>(
            materialName.c_str(), 
            Material(abs, nuFiss, fiss, chi, scatTable)));

	}
}

void MaterialLib::assignID(int id, const char* name){
    m_ids.insert(std::pair<int, const char*>(id, name));
    return;
}

};
