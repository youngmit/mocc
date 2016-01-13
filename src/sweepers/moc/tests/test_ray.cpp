#undef __STRICT_ANSI__
#undef _REENT_ONLY
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <cassert>
#include <string>

#include "pugixml.hpp"

#include "moc/ray_data.hpp"

#include "core/angular_quadrature.hpp"
#include "core/constants.hpp"
#include "core/core_mesh.hpp"
#include "core/global_config.hpp"

using namespace mocc;
using namespace mocc::moc;

using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE( testsimple )
{
    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file( "6x5.xml" );

    BOOST_CHECK( result );

    mocc::CoreMesh mesh( geom_xml );
    {

        // Test a few rays that starts on a corner, ends on a corner and crosses
        // a bunch of corners
        {
            Ray ray( Point2(0.0,1.0), Point2(4.0,5.0), 0, 0, 0, mesh );

            BOOST_CHECK_EQUAL( ray.cm_surf_fw(), 37 );
            BOOST_CHECK_EQUAL( ray.cm_cell_fw(), 6 );
            BOOST_CHECK_EQUAL( ray.cm_surf_bw(), 88 );
            BOOST_CHECK_EQUAL( ray.cm_cell_bw(), 27 );

            BOOST_CHECK_EQUAL( ray.nseg(), 12 );
            BOOST_CHECK_EQUAL( ray.ncseg(), 8 );

            // all of the segment lengths should be the same. I'm not testing
            // this too much in the general sense, since the tests for the pin
            // meshes should find most of these types of issues.
            real_t t = 1.0/3.0*sqrt(2);
            for( auto v: ray.seg_len() ) {
                BOOST_CHECK_CLOSE( v, t, 0.00001 );
            }

            std::vector<Surface> fw_surf =
            {
                Surface::EAST,
                Surface::NORTH,
                Surface::EAST,
                Surface::NORTH,
                Surface::EAST,
                Surface::NORTH,
                Surface::EAST,
                Surface::NORTH
            };

            std::vector<Surface> bw_surf =
            {
                Surface::WEST,
                Surface::SOUTH,
                Surface::WEST,
                Surface::SOUTH,
                Surface::WEST,
                Surface::SOUTH,
                Surface::SOUTH,
                Surface::WEST
            };

            VecI nseg = { 3, 0, 3, 0, 3, 0, 3, 0 };
            for( size_t i=0; i<ray.ncseg(); i++ ) {
                auto rcd = ray.cm_data()[i];
                BOOST_CHECK_EQUAL( rcd.fw, fw_surf[i] );
                BOOST_CHECK_EQUAL( rcd.bw, bw_surf[i] );
                BOOST_CHECK_EQUAL( rcd.nseg_fw, nseg[i] );
                BOOST_CHECK_EQUAL( rcd.nseg_bw, nseg[i] );
            }
        }

        {
            Ray ray( Point2(4.0,0.0), Point2(6.0,2.0), 0, 0, 0, mesh );

            BOOST_CHECK_EQUAL( ray.cm_surf_fw(), 89 );
            BOOST_CHECK_EQUAL( ray.cm_cell_fw(), 4 );
            BOOST_CHECK_EQUAL( ray.cm_surf_bw(), 43 );
            BOOST_CHECK_EQUAL( ray.cm_cell_bw(), 11 );

            BOOST_CHECK_EQUAL( ray.nseg(), 6 );
            BOOST_CHECK_EQUAL( ray.ncseg(), 4 );

        }

        {
            Ray ray( Point2(2.0,0.0), Point2(0.0,2.0), 0, 0, 0, mesh );
        }

        {
            Ray ray( Point2(6.0,3.0), Point2(4.0,5.0), 0, 0, 0, mesh );
        }

        {
            Ray ray( Point2(0.0,0.5), Point2(6.0,3.25), 0, 0, 0, mesh );
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
    }
    {
        Ray ray( Point2(3.78, 2.52), Point2(2.52, 3.78), 0, 0, 0, mesh );
    }
}
