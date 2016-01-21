#include "UnitTest++/UnitTest++.h"

#include "core/assembly.hpp"
#include "core/lattice.hpp"


TEST( assembly )
{
    std::string test_xml =
        "<mesh type=\"rec\" pitch=\"1.2\">"
        "<sub_x>1</sub_x>"
        "<sub_y>1</sub_y>"
        "</mesh>"
        ""
        "<material_lib path=\"c5g7.xsl\">"
        "   <material id=\"1\" name=\"UO2-3.3\" />"
        "</material_lib>"
        ""
        "<pin id=\"1\" mesh=\"1\">"
        "   1"
        "</pin>"
        ""
        "<lattice id=\"1\" nx=\"3\" ny=\"5\">"
        "   1 1 1"
        "   1 1 1"
        "   1 1 1"
        "   1 1 1"
        "   1 1 1"
        "</lattice>"
        ""
        "<assembly id=\"1\" np=\"5\" hz=\"3.14\">"
        "   <lattices>"
        "       1"
        "       1"
        "       1"
        "       1"
        "       1"
        "   </lattices>"
        "</assembly>"
        ""
        "<assembly id=\"1\" np=\"5\" hz=\"3.14\">"
        "   <lattices>"
        "       1"
        "       1"
        "       1"
        "       1"
        "       1"
        "   </lattices>"
        "</assembly>"
        "";

    pugi::xml_document xml;
    auto result = xml.load_string( test_xml.c_str() );

    if( !result ){
        return 1;
    }
}

int main(int, const char*[]) {
    return UnitTest::RunAllTests();
}
