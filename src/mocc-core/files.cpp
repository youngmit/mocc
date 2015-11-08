#include "files.hpp"

#include <iostream>
#include <string>
#include <string.h>

std::fstream LogFile;
std::fstream OutFile;

// A utility function for stripping the extension from the end of the command
// line argument and replacing with '.log'
void StartLogFile(const char* arg){
    std::string fname = arg;
	size_t pos = fname.rfind(".");
    std::string logname = fname.substr(0, pos);
    logname.append(".log");

	std::cout << "Logging output to: " << logname << std::endl << std::endl;
	LogFile.open(logname, std::fstream::out);
}

void StopLogFile(){
    LogFile.close();
}
