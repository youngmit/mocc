#include <string>
#include <cassert>

#include "pugixml.hpp"

#include "ray_data.hpp"
#include "angular_quadrature.hpp"
#include "core_mesh.hpp"
#include "constants.hpp"
#include "global_config.hpp"

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

    mocc::AngularQuadrature ang_quad( angquad_xml.child("ang_quad") );

    pugi::xml_document ray_xml;
    ray_xml.load_string("<rays spacing=\"0.01\" />");

    mocc::RayData ray_data( ray_xml.child("rays"), ang_quad, mesh );



    for( auto &plane_rays: ray_data ) {
        int iang = 0;
        double wsum = 0.0;
        mocc::VecF vol(mesh.n_reg(), 0.0);
        for( auto &angle_rays: plane_rays ) {

            wsum += ang_quad[iang].weight*2.0*PI;
            mocc::float_t space = ray_data.spacing(iang);
            mocc::float_t wt_ang = space * ang_quad[iang].weight * 2.0*PI;
            for( auto &ray: angle_rays ) {
                for( unsigned int iseg=0; iseg<ray.nseg(); iseg++ ) {
                    int ireg = ray.seg_index(iseg);
                    vol[ireg] += ray.seg_len(iseg) * wt_ang;
                }
            }
            
            iang++;
        }
        for( auto &v: vol ) {
            std::cout << v/0.1764 << std::endl;
        }
        std::cout << "angle weight sum: " << wsum << std::endl;
    }

    
    return 0;
}
