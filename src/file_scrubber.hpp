#pragma once
#include <string>
#include <fstream>

// A simple class that extends ifstream to provide a safe/easy way to scrub
// comments from an input file as it is read

class FileScrubber{
public:
	// Ininialize from a file name. Passees fName into the constructor of the 
	// underlying ifstream object.
	FileScrubber(const char* fName, const char* commentFlag);
	// Destructor doesnt really do anything right now.
	~FileScrubber();
	// Extend getline to only return lines which are not empty after removing
	// comments, stripped of comments.
	std::string getline();
	// Wrap the eof() function on the underlying ifstream object
	bool eof(){
		return m_stream.eof();
	}
	
private:
	std::ifstream m_stream;
	std::string m_flag;
};