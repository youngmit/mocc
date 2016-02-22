#include "files.hpp"

#include <iostream>
#include <string>
#include <string.h>

std::fstream LogFile;
std::fstream OutFile;

teestream LogScreen(std::cout, std::cout);

// A utility function for stripping the extension from the end of the command
// line argument and replacing with '.log'
void StartLogFile(std::string arg) {
    std::string logname = arg;
    logname.append(".log");

    std::cout << "Logging output to: " << logname << std::endl << std::endl;
    LogFile.open(logname, std::fstream::out);

    LogScreen = teestream(std::cout, LogFile);

    return;
}

void StopLogFile() {
    LogFile.close();
}
