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

#pragma once

#include <algorithm>
#include <string>
#include <vector>

// From http://stackoverflow.com/a/25829233

// trim from left
inline std::string &ltrim(std::string &s, const char *t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from right
inline std::string &rtrim(std::string &s, const char *t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from left & right
inline std::string &trim(std::string &s, const char *t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

// copying versions

inline std::string ltrim_copy(std::string s, const char *t = " \t\n\r\f\v")
{
    return ltrim(s, t);
}

inline std::string rtrim_copy(std::string s, const char *t = " \t\n\r\f\v")
{
    return rtrim(s, t);
}

inline std::string trim_copy(std::string s, const char *t = " \t\n\r\f\v")
{
    return trim(s, t);
}

/**
 * \brief Return a string, representing the ranges of a std::vector<bool> that
 * are true.
 */
std::string print_range(const std::vector<bool> &input);

// Sanitize a string: remove whitespace and cast to lowercase.
inline std::string &sanitize(std::string &s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return trim(s);
}

/**
 * \brief Given a string, nominally containing whitespace-delimited integers,
 * return a vector of those integers
 *
 * This will fail if the string contains characters that are not numererals or
 * whitespace.
 */
template <typename T> std::vector<T> explode_string(std::string data);
