#include "geometry_output.hpp"

#include <fstream>
#include <string>

#include "pugixml.hpp"

#include "error.hpp"

using std::ofstream;
using std::endl;
using std::string;

namespace mocc { namespace aux {

    void output_geometry( const pugi::xml_node &input, const CoreMesh &mesh )
    {
        if( input.empty() ){
            throw EXCEPT("No input for geometry output.");
        }

        if( input.attribute("file").empty() ) {
            throw EXCEPT("No \"file\" attribute specified.");
        }
        string file = input.attribute("file").value();

        int plane = input.attribute("plane").as_int(0);
        if( (plane < 0) || (plane >= (int)mesh.nz() ) ) {
            throw EXCEPT("Invalid plane specified.");
        }

        ofstream out(file);

        // boilerplate
        out << "import cairo as cr" << endl;
        out << "import math" << endl;
        out << "import rays" << endl;
        out << endl;
        out << "twopi = math.pi*2" << endl;
        out << endl;
        out << "# set this to whichever angle of ray you want to show." 
                " Negative value to disable." << endl;
        out << "angle = -1" << endl;
        out << endl;
        out << "mesh_lines = []" << endl;
        out << endl;
        out << "core_dims = [" << mesh.hx_core() << ", " 
                               << mesh.hy_core() << "]" << endl;
        out << "" << endl;
        out << "surface = cr.PDFSurface(\"geometry.pdf\", 720, 720)" << endl;
        out << "ctx = cr.Context(surface)" << endl;
        out << "ctx.scale(720/core_dims[0], -720/core_dims[1])" << endl;
        out << "ctx.translate(0, -core_dims[1])" << endl;
        out << "" << endl;


        
        // do all of the mesh stuff
        out << "ctx.set_line_width(0.001)" << endl;
        out << "" << endl;
        out << "ctx.set_source_rgb(0, 0, 0)" << endl;
        out << "" << endl;

        for( auto l: mesh.lines() ) {
            out << "mesh_lines.append(" << l << ")" << endl;
        }

        // Draw the core lines
        out << "" << endl;
        out << "for l in mesh_lines:" << endl;
        out << "    p1 = l[0]" << endl;
        out << "    p2 = l[1]" << endl;
        out << "    ctx.move_to(p1[0], p1[1])" << endl;
        out << "    ctx.line_to(p2[0], p2[1])" << endl;

        out << endl;

        // Draw the pin meshes
        int ipin = 0;
        for( auto pin=mesh.begin( plane );
                pin!=mesh.end( plane ); ++pin )
        {
            out << "print \"drawing pin \" + str(" << ipin << ")" << endl;
            const PinMesh &pm = (*pin)->mesh();
            Point2 origin = mesh.pin_origin(ipin);

            out << "ctx.translate(" << origin.x << ", " 
                                    << origin.y << ")" << endl;

            out << pm.draw() << endl;
            out << "ctx.translate(" << -origin.x << ", " 
                                    << -origin.y << ")" << endl << endl;
            ipin++;
        }

        out << endl;

        // Do ray output
        out << "if angle >= 0:" << endl;
        out << "    ctx.set_source_rgb(0, 0, 1)" << endl;
        out << "    rays.draw_rays(ctx, angle)" << endl;
        out << "" << endl;

        out << "surface.finish()" << endl;
        out << "" << endl;

        return;
    }

} } // namespaces
