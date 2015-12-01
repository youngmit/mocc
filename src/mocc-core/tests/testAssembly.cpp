#undef __STRICT_ANSI__
#undef _REENT_ONLY
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "mocc-core/assembly.hpp"
#include "mocc-core/lattice.hpp"


BOOST_AUTO_TEST_CASE( testall )
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
    auto result = xml.load_string( xml );

    if( !result ){
        return 1;
    }


}
