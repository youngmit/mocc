// Global definitions for output files. This keeps us from having to pass
// instances of the log and output files all over the place.
#pragma once

#include <fstream>

extern std::fstream LogFile;
extern std::fstream OutFile;
extern std::string CaseName;

void StartLogFile(const char* arg);

void StopLogFile();
