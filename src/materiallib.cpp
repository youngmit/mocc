#include "materiallib.hpp"
#include "material.hpp"
#include "error.hpp"
#include <string>
#include <sstream>
#include <regex>
#include <iostream>

using std::stringstream;
using std::string;

namespace mocc{

MaterialLib::MaterialLib(FileScrubber &input){
	string line;
	// Read the first three lines to extract the library header
	// Description
	m_description = input.getline();
	// Number of groups and number of materials
	{
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
		
		std::cout << "0" << std::endl;
		// std::regex headExp("^\\s*XSMACRO\\s+([^\\s]+)\\s+([0-9]+)\\s*$");
		try{
			std::regex headExp("\\s XSMACRO", std::regex_constants::extended);
			std::cout << "a" << std::endl;
			std::smatch results;
			std::cout << "b" << std::endl;
			regex_match(line, results, headExp);
			std::cout << "c" << std::endl;
			std::cout << results[1].str() << std::endl;
		} catch(std::regex_error e) {
			std::cout << "regex error code: ";
			switch(e.code()){
				case std::regex_constants::error_collate:
					std::cout << "error_collate";
					break;
				case std::regex_constants::error_ctype:
					std::cout << "error_ctype";
					break;
				case std::regex_constants::error_escape:
					std::cout << "error_escape";
					break;
				case std::regex_constants::error_backref:
					std::cout << "error_backref";
					break;
				case std::regex_constants::error_brack:
					std::cout << "error_brack";
					break;
				case std::regex_constants::error_paren:
					std::cout << "error_paren";
					break;
				case std::regex_constants::error_brace:
					std::cout << "error_brace";
					break;
				case std::regex_constants::error_badbrace:
					std::cout << "error_badbrace";
					break;
				default:
					break;
			}
			std::cout << std::endl;
		}
		
		
	}
}

};