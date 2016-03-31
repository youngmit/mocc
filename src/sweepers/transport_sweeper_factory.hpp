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

#include "pugixml.hpp"

#include "transport_sweeper.hpp"

namespace mocc {
    /**
    * Peek inside a \<sweeper\> tag to look at the \c type attribute, then
    * generate a \ref TransportSweeper of the appropriate type using the passed
    * XML node and \ref CoreMesh.
    */
    UP_Sweeper_t TransportSweeperFactory( const pugi::xml_node &input,
        const CoreMesh& mesh );
}
