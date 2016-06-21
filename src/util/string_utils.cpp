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
#include <iostream>
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
        } else {
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

    if (data.size() == 0) {
        return std::vector<T>();
    }

    std::vector<T> out;

    // Store the string as a stream and try to read all entries
    std::stringstream inBuf(data);
    T i;
    while (!inBuf.eof()) {
        inBuf >> i;
        if (!inBuf.fail()) {
            out.push_back(i);
        } else {
            std::stringstream msg;
            msg << "Malformed data: " << data;
            throw EXCEPT(msg.str());
        }
    }

    // Make sure everything went okay
    if (inBuf.fail()) {
        throw EXCEPT("Trouble reading data");
    }

    return out;
}

/**
 * \brief Break a string with matched curly braces ({ }).
 *
 * Return a vector of strings containing the substrings enclosed in each pair of
 * braces. Any value
 *
 * Throws an exception under the following circumstances:
 *  - If there are any characters that are not whitespace, numerals, or braces
 *  - If there are non-whitespace characters that are not enclosed in a brace
 * pair
 *  - If there are un-matched braces
 *  - If the brace depth exceeds one
 *
 */
std::vector<std::vector<int>> explode_braces(std::string data)
{
    auto is_invalid_char = [](char c) {
        return !(isspace(c) || isdigit(c) || (c == '{') || (c == '}'));
    };

    // first, make sure there are no non-[numerals|whitespace|braces]
    bool good =
        std::find_if(data.begin(), data.end(), is_invalid_char) == data.end();
    if (!good) {
        throw EXCEPT("Malformed data");
    }

    // Make sure that all of the braces match.
    int brace_depth = 0;
    for (const auto &c : data) {
        if (c == '{') {
            brace_depth++;
            if (brace_depth > 1) {
                break;
            }
            continue;
        }
        if (c == '}') {
            brace_depth--;
            continue;
        }
    }
    if (brace_depth != 0) {
        throw EXCEPT("Brace mismatch");
    }

    // Read the actual data
    std::vector<std::vector<int>> out;
    size_t pos_stt = 0;
    size_t pos_stp = 0;
    size_t pos     = 0;
    while (true) {
        pos_stt = data.find('{', pos);
        pos_stp = data.find('}', pos_stt);
        if (pos_stt - pos > 0) {
            std::vector<int> naked_ints =
                explode_string<int>(data.substr(pos, pos_stt - pos - 1));
            if (naked_ints.size() > 0) {
                for (const auto &v : naked_ints) {
                    out.push_back(std::vector<int>(1, v));
                }
            }
        }
        if (pos_stt != std::string::npos) {
            out.push_back(explode_string<int>(
                data.substr(pos_stt + 1, pos_stp - pos_stt - 1)));
            if (out.back().size() == 0) {
                throw EXCEPT("Empty brackets");
            }
        } else {
            break;
        }
        pos = pos_stp + 1;
    }

    return out;
}

template std::vector<int> explode_string(std::string data);
template std::vector<mocc::real_t> explode_string(std::string data);
