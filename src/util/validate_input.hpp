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

#include <string>
#include <vector>
#include "util/error.hpp"
#include "util/pugifwd.hpp"

namespace mocc {
/**
 * \brief Loop over all attributes of a passed XML node and check validity
 * 
 * \param input the XML node to validate
 * \param recognized_attributes a vector containing the names of the attributes
 * that are considered valid for the passed XML node.
 *
 * This will generate a warning if any of the attributes on the passed XML node
 * cannot be found in the vector of recognized attributes. This function does
 * not check validity of the attribute values themselves. It is assumed the the
 * client code, knowing which attributes to parse, performs its own checks of
 * the validity of attribute values.
 */
bool validate_input(const pugi::xml_node &input,
                    const std::vector<std::string> recognized_attributes);
}
