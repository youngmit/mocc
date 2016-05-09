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

#include <cassert>
#include <iostream>
#include <memory>
#include "pugixml.hpp"

/**
 * \brief Make an XML document inline, from a string containing the XML
 */
auto inline_xml( const char* input ) {
    std::unique_ptr<pugi::xml_document>
        doc(std::make_unique<pugi::xml_document>());
    pugi::xml_parse_result result = doc->load_string( input );

    if(!result) {
        std::cout << result.description() << std::endl;
    }
    assert( result );

    return doc;
}

/**
 * \brief Make an XML document inline, from a filename
 */
auto inline_xml_file( const char* input ) {
    std::unique_ptr<pugi::xml_document>
        doc(std::make_unique<pugi::xml_document>());
    pugi::xml_parse_result result = doc->load_file( input );

    if(!result) {
        std::cout << result.description() << std::endl;
    }
    assert( result );

    return doc;
}
