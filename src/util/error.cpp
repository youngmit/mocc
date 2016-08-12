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

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "error.hpp"
#include "files.hpp"

#include "util/stacktrace.hpp"

using std::string;

namespace mocc {

std::unordered_map<std::string, Warning> Warnings;

void Error(const char *msg)
{
    std::cerr << "ERROR: " << msg << std::endl;
    LogFile << "ERROR: " << msg << std::endl;
    exit(EXIT_FAILURE);
}

void Warn(const std::string &msg)
{
    auto it = Warnings.find(msg);
    if (it == Warnings.end()) {
        Warnings.emplace(msg, msg);
        LogScreen << "WARNING: " << msg << std::endl;
    }
    else {
        it->second.count++;
    }
}

void Fail(Exception e)
{
    std::cerr << e.what();
    LogFile << e.what();
    exit(EXIT_FAILURE);
}

Exception::Exception(Info info)
    : info_(info)
{
    std::stringstream ret;
    ret << info_.file << ":" << info.line << " in " << info.func << std::endl;
    ret << info_.msg << std::endl;
    print_message_ = ret.str();

    return;
}

Exception::Exception(Info info, const Exception &parent)
    : Exception(info)
{
    print_message_.append(parent.what());

    return;
}

const char *Exception::what() const noexcept
{
    return print_message_.c_str();
}

std::ostream &operator<<(std::ostream &os, const Warning &warn)
{
    os << warn.count << ": " << warn.description;
    return os;
}
}
