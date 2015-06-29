#include "file_scrubber.hpp"
#include "string_utils.hpp"
#include <stdio.h>
#include <iostream>

using std::ifstream;


FileScrubber::FileScrubber(const char* fName, const char* commentFlag):
	m_stream(fName),
	m_flag(commentFlag){}
	
FileScrubber::~FileScrubber(){}

std::string FileScrubber::getline(){
	while(!m_stream.eof()){
		std::string line;
		std::getline(m_stream, line);
		// Strip the comments
		int commentPos = line.find(m_flag, 0);
		if(commentPos != std::string::npos){
			line.erase(line.find(m_flag, 0), std::string::npos);
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