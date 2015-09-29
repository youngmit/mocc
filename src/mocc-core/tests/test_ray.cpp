#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <cassert>
#include <string>

#include "pugixml.hpp"

#include "angular_quadrature.hpp"
#include "constants.hpp"
#include "core_mesh.hpp"
#include "global_config.hpp"
#include "ray_data.hpp"

using namespace mocc;

using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE( testsimple )
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file( "6x5.xml" );

    BOOST_CHECK( result );
    
    mocc::CoreMesh mesh( geom_xml );
    {
        
        {
            Ray ray( Point2(0.0,1.0), Point2(4.0,5.0), 0, 0, 0, mesh );

            //BOOST_CHECK(ray.);
            cout << ray.cm_data()[3].fw << endl;
            
            cout << "size of cm data element: " << sizeof(ray.cm_data().front()) << endl;
            for( auto rcd: ray.cm_data() ) {
                cout << rcd << endl;
            }
cout << endl;
        }

        {
            Ray ray( Point2(4.0,0.0), Point2(6.0,2.0), 0, 0, 0, mesh );
        }

        {
            Ray ray( Point2(2.0,0.0), Point2(0.0,2.0), 0, 0, 0, mesh );
        }

        {
            Ray ray( Point2(6.0,3.0), Point2(4.0,5.0), 0, 0, 0, mesh );
        }

        {
            Ray ray( Point2(0.0,0.5), Point2(6.0,3.25), 0, 0, 0, mesh );
            for( auto rcd: ray.cm_data() ) {
                cout << rcd << endl;
            }
cout << endl;
        }


    }
}

BOOST_AUTO_TEST_CASE( testall )
{

    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file( "square.xml" );
    
    mocc::CoreMesh mesh( geom_xml );

    pugi::xml_document angquad_xml;
    result = 
        angquad_xml.load_string("<ang_quad type=\"ls\" order=\"4\" />");
    // Make a nasty ray to exercise the coarse indexing
    {
        Ray ray( Point2(1.26, 0.0), Point2(3.78, 2.52), 0, 0, 0, mesh );
    }
    {
        Ray ray( Point2(1.26, 0.0), Point2(0.0, 1.26), 0, 0, 0, mesh );
    }
    {
        Ray ray( Point2(0.0, 1.26), Point2(2.52, 3.78), 0, 0, 0, mesh );
            for( auto rcd: ray.cm_data() ) {
                cout << rcd << endl;
            }
    }
    {
        Ray ray( Point2(3.78, 2.52), Point2(2.52, 3.78), 0, 0, 0, mesh );
    }
}
