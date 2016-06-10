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

#include "string_utils.hpp"

#include <ctype.h>
#include <sstream>
#include "util/error.hpp"
#include "util/global_config.hpp"

using mocc::Exception;

std::string print_range(const std::vector<bool> &input)
{
    std::stringstream output;
    bool left      = false;
    int left_bound = 0;
    std::vector<std::pair<int, int>> ranges;
    for (int i = 0; i < (int)input.size(); i++) {
        if (input[i] && !left) {
            left       = true;
            left_bound = i;
        }
        if (left && !input[i]) {
            left = false;
            ranges.emplace_back(left_bound, i - 1);
        }
        if (left && input[i] && (i == (int)input.size() - 1)) {
            ranges.emplace_back(left_bound, i);
        }
    }

    for (auto r : ranges) {
        if (r.first != r.second) {
            output << r.first << "-" << r.second;
        }
        else {
            output << r.first;
        }

        if (r != ranges.back()) {
            output << ", ";
        }
    }
    return output.str();
}

template <typename T> std::vector<T> explode_string(std::string data)
{
    sanitize(data);
    auto is_invalid_char = [](char c) { return !(isspace(c) || isdigit(c)); };

    // first, make sure there are no non-[numerals|whitespace]
    bool good =
        std::find_if(data.begin(), data.end(), is_invalid_char) == data.end();
    if (!good) {
        throw EXCEPT("Malformed data");
    }

    std::vector<T> out;

    // Store the string as a stream and try to read all entries
    std::stringstream inBuf(data);
    T i;
    while (!inBuf.eof()) {
        inBuf >> i;
        out.push_back(i);
    }

    // Make sure everything went okay
    if (inBuf.fail()) {
        throw EXCEPT("Trouble reading data");
    }

    return out;
}

template std::vector<int> explode_string(std::string data);
template std::vector<mocc::real_t> explode_string(std::string data);
