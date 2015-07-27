#include <cassert>
#include "pugixml.hpp"

pugi::xml_document inline_xml( const char* input ) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string( input );

    assert( result );

    return doc;
}
