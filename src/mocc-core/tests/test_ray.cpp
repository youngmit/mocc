#include <string>
#include <cassert>

#include "pugixml.hpp"

#include "ray_data.hpp"
#include "angular_quadrature.hpp"
#include "core_mesh.hpp"
#include "constants.hpp"
#include "global_config.hpp"

using namespace mocc;

using std::cout;
using std::endl;

int main() {

    pugi::xml_document geom_xml;
    pugi::xml_parse_result result = geom_xml.load_file( "square.xml" );
    
    mocc::CoreMesh mesh( geom_xml );

    pugi::xml_document angquad_xml;
    result = 
        angquad_xml.load_string("<ang_quad type=\"ls\" order=\"4\" />");

    if ( !result ) {
        return 1;
    }


    // Make a nasty ray to exercise the coarse indexing
    {
        Ray ray( Point2(1.26, 0.0), Point2(3.78, 2.52), 0, 0, 0, mesh );
        assert( ray.cm_surf(0) == 21 );
        assert( ray.cm_surf(1) == 10 );
        assert( ray.cm_surf(2) == 11 );
        assert( ray.cm_surf(3) == 30 );
        assert( ray.cm_surf(4) == 16 );
    }
    {
        Ray ray( Point2(0.0, 1.26), Point2(1.26, 0.0), 0, 0, 0, mesh );
        assert( ray.cm_surf(0) == 21 );
        assert( ray.cm_surf(1) == 9 );
    }
    {
        Ray ray( Point2(0.0, 1.26), Point2(2.52, 3.78), 0, 0, 0, mesh );
        assert( ray.cm_surf(0) == 9 );
        assert( ray.cm_surf(1) == 22 );
        assert( ray.cm_surf(2) == 14 );
        assert( ray.cm_surf(3) == 27 );
        assert( ray.cm_surf(4) == 32 );
        assert( ray.cm_surf(5) == 19 );
    }
    {
        Ray ray( Point2(2.52, 3.78), Point2(3.78, 2.52), 0, 0, 0, mesh );
        assert( ray.cm_surf(0) == 16 );
        assert( ray.cm_surf(1) == 31 );
        assert( ray.cm_surf(2) == 32 );
    }

    
    return 0;
}
