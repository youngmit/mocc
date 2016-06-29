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

#include "validate_input.hpp"

#include <algorithm>
#include <sstream>
#include "pugixml.hpp"

namespace mocc {
bool validate_input(const pugi::xml_node &input,
                    const std::vector<std::string> recognized_attributes)
{
    bool good = true;

    for (const auto attrib : input.attributes()) {
        if (std::find(recognized_attributes.begin(),
                      recognized_attributes.end(),
                      attrib.name()) == recognized_attributes.end()) {
            good = false;
            std::stringstream msg;
            msg << "Unrecognized attribute defined on <" << input.name()
                << "> tag: " << attrib.name();
            Warn(msg.str());
        }
    }

    return good;
}
}
