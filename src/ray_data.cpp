#include "ray_data.hpp"

#include <iostream>
#include <cmath>

#include "error.hpp"

using std::cout;
using std::endl;

namespace mocc {

    // Construct a ray from a starting point and angle
    Ray::Ray( Point2 p, Angle ang ) {

        return;
    }


    RayData::RayData( const pugi::xml_node &input, 
            const AngularQuadrature &ang_quad,
            const CoreMesh &mesh ):
    ang_quad_(ang_quad)
    {
cout << "Generating Ray Data" << endl;
        // Make sure we have reasonable input
        if ( input.empty() ) {
            Error("No input privided for ray spacing.");
        }

        // Get the optimal ray spacing
        float_t opt_spacing = input.attribute("spacing").as_float(-1.0);
        if( opt_spacing <= 0.0 ) {
            Error("Failed to read valid ray spacing.");
        }

        // Figure out modular angles and spacings
        float_t hx = mesh.hx();
        float_t hy = mesh.hy();

        for (auto ang_it = ang_quad_.octant(1); 
                ang_it != ang_quad_.octant(3); ++ang_it) {
            Angle ang = *ang_it;

            
            int Nx = ceil(hx/opt_spacing*fabs(sin(ang.alpha)));
            int Ny = ceil(hy/opt_spacing*fabs(cos(ang.alpha)));
            Nx += Nx%2;
            Ny += Ny%2; 


            float_t new_alpha = atan(hy*Nx/(hx*Ny));
            float_t actual_spacing = 0;
cout << "Number of rays: " << Nx << " " << Ny << endl;
cout << "Old/New alpha: " << ang.alpha << "/" << new_alpha << endl;
            
        }

        //
        
        // Trace rays

        // TODO: adjust ray lengths to correct FSR volume
    }


}
