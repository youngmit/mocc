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

#include "file_scrubber.hpp"

#include <stdio.h>
#include <iostream>

#include "string_utils.hpp"
#include "error.hpp"

using std::ifstream;

namespace mocc {
    FileScrubber::FileScrubber(const char* fName, const char* commentFlag):
        stream_(fName),
        flag_(commentFlag)
    {
        if( !stream_.good() ) {
            throw EXCEPT("Failed to open file");
        }
    }

    std::string FileScrubber::getline() {
        while(!stream_.eof()) {
            std::string line;
            std::getline(stream_, line);
            // Strip the comments
            size_t commentPos = line.find(flag_, 0);
            if(commentPos != std::string::npos) {
                line.erase(commentPos, std::string::npos);
            }
            // Remove whitespace
            line = trim(line);
            // If the result isnt empty, return the line.
            if(!line.empty()) {
                return line;
            }
        }
        return "";
    }
}
