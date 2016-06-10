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

// This exists as the entry point into all pin mesh types, also providing a
// factory method for generating deferred-type pin mesh objects

#include "pin_mesh.hpp"

#include <iostream>
#include <sstream>

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"

namespace mocc {
// Determine which type of pin to create from an XML object, produce a mesh
// of the appropriate type and return a shared pointer to the object.
PinMesh *PinMeshFactory(const pugi::xml_node &input)
{
    PinMesh *pm = nullptr;

    // Extract the type of mesh to make
    std::string type = input.attribute("type").value();

    if (type == "cyl") {
        pm = new PinMesh_Cyl(input);
    }
    else if (type == "rect") {
        pm = new PinMesh_Rect(input);
    }
    else {
        // I don't recognize the mesh type, error out.
        std::stringstream err;
        err << "Unrecognized mesh type for mesh ID: "
            << input.attribute("id").value();
        throw EXCEPT(err.str());
    }

    return pm;
}

std::map<int, UP_PinMesh_t> ParsePinMeshes(const pugi::xml_node &input)
{
    std::map<int, UP_PinMesh_t> pin_meshes;
    for (pugi::xml_node mesh = input.child("mesh"); mesh;
         mesh                = mesh.next_sibling("mesh")) {
        LogFile << "Parsing new pin mesh: ID=" << mesh.attribute("id").value()
                << std::endl;
        UP_PinMesh_t pm(PinMeshFactory(mesh));
        int id = pm->id();
        if (pin_meshes.find(id) != pin_meshes.end()) {
            std::stringstream msg;
            msg << "Duplicate pin mesh ID (" << id << ")";
            throw EXCEPT(msg.str());
        }
        pin_meshes.emplace(id, std::move(pm));
    }

    return pin_meshes;
}
}
