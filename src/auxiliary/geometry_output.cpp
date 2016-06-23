/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "geometry_output.hpp"
#include <fstream>
#include <string>
#include "pugixml.hpp"
#include "util/error.hpp"

using std::ofstream;
using std::string;

namespace mocc {
namespace aux {

void output_geometry(const pugi::xml_node &input, const CoreMesh &mesh)
{
    if (input.empty()) {
        throw EXCEPT("No input for geometry output.");
    }

    string file;
    if (input.attribute("file").empty()) {
        file = "geom.py";
    }
    else {
        file = input.attribute("file").value();
    }

    int plane = input.attribute("plane").as_int(0);
    if ((plane < 0) || (plane >= (int)mesh.nz())) {
        throw EXCEPT("Invalid plane specified.");
    }

    ofstream out(file);

    // boilerplate
    out << "import cairo as cr" << std::endl;
    out << "import math" << std::endl;
    out << "import rays" << std::endl;
    out << std::endl;
    out << "twopi = math.pi*2" << std::endl;
    out << std::endl;
    out << "# set this to whichever angle of ray you want to show."
           " Negative value to disable."
        << std::endl;
    out << "angle = -1" << std::endl;
    out << std::endl;
    out << "mesh_lines = []" << std::endl;
    out << std::endl;
    out << "core_dims = [" << mesh.hx_core() << ", " << mesh.hy_core() << "]"
        << std::endl;
    out << "" << std::endl;
    out << "surface = cr.PDFSurface(\"geometry.pdf\", 720, 720)" << std::endl;
    out << "ctx = cr.Context(surface)" << std::endl;
    out << "ctx.scale(720/core_dims[0], -720/core_dims[1])" << std::endl;
    out << "ctx.translate(0, -core_dims[1])" << std::endl;
    out << "" << std::endl;

    // do all of the mesh stuff
    out << "ctx.set_line_width(0.001)" << std::endl;
    out << "" << std::endl;
    out << "ctx.set_source_rgb(0, 0, 0)" << std::endl;
    out << "" << std::endl;

    for (auto l : mesh.lines()) {
        out << "mesh_lines.append(" << l << ")" << std::endl;
    }

    // Draw the core lines
    out << "" << std::endl;
    out << "for l in mesh_lines:" << std::endl;
    out << "    p1 = l[0]" << std::endl;
    out << "    p2 = l[1]" << std::endl;
    out << "    ctx.move_to(p1[0], p1[1])" << std::endl;
    out << "    ctx.line_to(p2[0], p2[1])" << std::endl;

    out << std::endl;

    // Draw the pin meshes
    int ipin = 0;
    for (auto pin = mesh.begin(plane); pin != mesh.end(plane); ++pin) {
        out << "print \"drawing pin \" + str(" << ipin << ")" << std::endl;
        const PinMesh &pm = (*pin)->mesh();
        Point2 origin     = mesh.pin_origin(ipin);

        out << "ctx.translate(" << origin.x << ", " << origin.y << ")" << std::endl;

        out << pm.draw() << std::endl;
        out << "ctx.translate(" << -origin.x << ", " << -origin.y << ")" << std::endl
            << std::endl;
        ipin++;
    }

    out << std::endl;

    // Do ray output
    out << "if angle >= 0:" << std::endl;
    out << "    ctx.set_source_rgb(0, 0, 1)" << std::endl;
    out << "    rays.draw_rays(ctx, angle)" << std::endl;
    out << "" << std::endl;

    out << "surface.finish()" << std::endl;
    out << "" << std::endl;

    return;
}
}
} // namespaces
