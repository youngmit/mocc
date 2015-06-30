// Global definitions for output files. This keeps us from having to pass
// instances of the log and output files all over the place.
#pragma once

#include <fstream>

extern std::fstream LogFile;
extern std::fstream OutFile;

// A utility function for stripping '.xml' from the end of the command line
// argument and replacing with '.log'
void StartLogFile(const char* arg);

void StopLogFile();
