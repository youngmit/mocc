#include "UnitTest++/UnitTest++.h"

#include "pugixml.hpp"

#include "xs_mesh_homogenized.hpp"

using namespace mocc;

TEST( xsmeshhom ) {
    {
        pugi::xml_document geom_xml;
        pugi::xml_parse_result result = geom_xml.load_file( "square.xml" );
        CHECK( result );

        CoreMesh mesh( geom_xml );

        XSMeshHomogenized xs_mesh( mesh );

        H5Node h5f( "xsmesh_1.h5", H5Access::WRITE );
        auto g = h5f.create_group("xs_mesh");
        xs_mesh.output(g.get());
    }
    {
        pugi::xml_document geom_xml;
        pugi::xml_parse_result result = geom_xml.load_file( "square_2.xml" );
        CHECK( result );

        CoreMesh mesh( geom_xml );

        XSMeshHomogenized xs_mesh( mesh );

        H5Node h5f( "xsmesh_2.h5", H5Access::WRITE );
        xs_mesh.output(h5f.get());
    }
    
}

TEST( fromdata ) {
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file( "stack.xml" );
    CHECK( result );
    if( !result ) {
        std::cout << result.description() << std::endl;
    }

    CoreMesh mesh( geom_xml );

    pugi::xml_document xsmesh_xml;
    std::string xml = 
        "<data file=\"xsmesh_1.h5\" top_plane=\"3\">"
        "<data file=\"xsmesh_2.h5\" top_plane=\"11\">"
        "";
    xsmesh_xml.load_string( xml.c_str() );

    XSMeshHomogenized xs_mesh( mesh, xsmesh_xml );

    
}

int main() {
    return UnitTest::RunAllTests();
}
