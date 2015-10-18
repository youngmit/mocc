#include "files.hpp"

#include <iostream>
#include <string>
#include <string.h>

std::fstream LogFile;
std::fstream OutFile;

// A utility function for stripping '.xml' from the end of the command line
// argument and replacing with '.log'
void StartLogFile(const char* arg){
	// See if our filename conforms to standard *.xml format.
	int len = strlen(arg);
	int pos = len;
	char* logname = new char[len+4];
	if (len > 4){
		const char* ext = &arg[len-4];
		if (strncmp(ext, ".xml", 4) == 0){
			pos = len-4;
		}
	}
	strncpy(logname, arg, len);
	
	// Actually add the .log
	strcpy(&logname[pos], ".log");
	std::cout << "Logging output to: " << logname << std::endl << std::endl;
	LogFile.open(logname, std::fstream::out);

    delete[] logname;
}

void StopLogFile(){
    LogFile.close();
}
