/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "files.hpp"

#include <iostream>
#include <string>
#include <string.h>

std::fstream LogFile;
std::fstream OutFile;

onullstream NullStream;

TeeStream LogScreen(std::cout, NullStream);

// A utility function for stripping the extension from the end of the command
// line argument and replacing with '.log'
void StartLogFile(std::string arg) {
    std::string logname = arg;
    logname.append(".log");

    std::cout << "Logging output to: " << logname << std::endl << std::endl;
    LogFile.open(logname, std::fstream::out);

    LogScreen.reset(std::cout, LogFile);

    return;
}

void StopLogFile() {
    LogFile.close();
}
