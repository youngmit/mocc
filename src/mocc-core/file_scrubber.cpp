#include "file_scrubber.hpp"

#include <stdio.h>
#include <iostream>

#include "string_utils.hpp"
#include "error.hpp"

using std::ifstream;

namespace mocc {
    FileScrubber::FileScrubber(const char* fName, const char* commentFlag):
    	stream_(fName),
    	flag_(commentFlag){
       
            if( !stream_.good() ) {
                Error("Failed to open file.");
            }
        
        }
    	
    FileScrubber::~FileScrubber(){}
    
    std::string FileScrubber::getline(){
    	while(!stream_.eof()){
    		std::string line;
    		std::getline(stream_, line);
    		// Strip the comments
    		unsigned int commentPos = line.find(flag_, 0);
    		if(commentPos != std::string::npos){
    			line.erase(commentPos, std::string::npos);
    		}
    		// Remove whitespace
    		line = trim(line);
    		// If the result isnt empty, return the line.
    		if(!line.empty()){
    			return line;
    		}
    	}
    	return "";
    }
}
