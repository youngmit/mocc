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

#include "assembly.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include "pugixml.hpp"

#include "error.hpp"
#include "string_utils.hpp"

using std::string;
using std::endl;
using std::cout;

namespace mocc {
Assembly::Assembly(const pugi::xml_node &input,
                   const std::map<int, UP_Lattice_t> &lattices)
{
    // Parse assembly ID
    id_ = input.attribute("id").as_int(0);
    if (id_ == 0) {
        throw EXCEPT("Invalid assembly ID.");
    }

    // Parse number of planes
    nz_ = input.attribute("np").as_int(0);
    if (nz_ == 0) {
        throw EXCEPT("Invalid number of planes (nz) when parsing "
                     "assembly.");
    }

    // Parse plane heights (scalar form)
    bool scalar_hz = false;
    real_t hz      = input.attribute("hz").as_float(0.0f);
    if (hz > 0.0f) {
        scalar_hz = true;
        // Fill the hz vector with all entries the same.
        dz_ = VecF(nz_, hz);
    }

    // Parse plane heights (array form)
    if (auto hz_in = input.child("hz")) {
        if (scalar_hz) {
            // hz is over-defined
            throw EXCEPT("Plane heights are over-specified for assembly.");
        }
        string hzs = hz_in.child_value();
        dz_        = explode_string<real_t>(hzs);

        // Lattice dimensions are read as top-to-bottom, but stored as
        // bottom-to-top.
        std::reverse(dz_.begin(), dz_.end());
    }

    // Parse lattice IDs
    if (input.child("lattices").empty()) {
        throw EXCEPT("No lattices specified!");
    }
    {
        string lat_str   = input.child("lattices").child_value();
        auto lattice_ids = explode_string<int>(lat_str);

        for (int lat_id : lattice_ids) {
            if (lattices.count(lat_id) > 0) {
                lattices_.push_back(lattices.at(lat_id).get());
            }
            else {
                throw EXCEPT("Unrecognized lattice ID in assembly.");
            }
        }
        if (lattices_.size() != nz_) {
            throw EXCEPT("Incorrect number of lattices specified for "
                         "assembly.");
        }
    }
    // Lattice IDs are read from the top down, but are stored from the
    // bottom up. Reverse the vector we just populated.
    std::reverse(lattices_.begin(), lattices_.end());

    // Store lattice dimensions
    hx_ = lattices_[0]->hx();
    hy_ = lattices_[0]->hy();

    // Make sure that all of the lattices are the same size
    for (const auto &lat : lattices_) {
        if (!lat->compatible(*lattices_[0])) {
            throw EXCEPT("Lattices in Assembly are not compatible.");
        }
    }

    // Store the total numver of FSRs and XS regions in the assembly
    n_reg_   = 0;
    n_xsreg_ = 0;
    for (auto &l : lattices_) {
        n_reg_ += l->n_reg();
        n_xsreg_ += l->n_xsreg();
    }

    return;
}

Assembly::~Assembly()
{
    return;
}

std::map<int, UP_Assembly_t>
ParseAssemblies(const pugi::xml_node &input,
                const std::map<int, UP_Lattice_t> lattices)
{
    std::map<int, UP_Assembly_t> assemblies;

    for (pugi::xml_node asy = input.child("assembly"); asy;
         asy                = asy.next_sibling("assembly")) {
        UP_Assembly_t asy_p(new Assembly(asy, lattices));
        int id = asy_p->id();
        if (assemblies.find(id) != assemblies.end()) {
            std::stringstream msg;
            msg << "Duplicate assembly ID (" << id << ")";
            throw EXCEPT(msg.str());
        }
        assemblies.emplace(asy_p->id(), std::move(asy_p));
    }

    return assemblies;
}
}
