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

#pragma once

#include <map>
#include <memory>
#include <vector>

#include "global_config.hpp"
#include "material_lib.hpp"
#include "pin_mesh.hpp"

namespace mocc {
/**
 * The Pin class is a concrete instantiaion of a physical pin. It
 * essentially applies materials to regions of a PinMesh. Nothin' fancy.
 */
class Pin {
public:
    Pin(const pugi::xml_node &input, const std::map<int, UP_PinMesh_t> &meshes,
        const MaterialLib &mat_lib);
    ~Pin()
    {
        return;
    }

    const PinMesh &mesh() const
    {
        return *pin_mesh_;
    }

    int id() const
    {
        return id_;
    }

    int n_reg() const
    {
        return pin_mesh_->n_reg();
    }

    int mesh_id() const
    {
        return pin_mesh_->id();
    }

    real_t vol() const
    {
        return pin_mesh_->vol();
    }

    const VecF &vols() const
    {
        return pin_mesh_->vols();
    }

    const VecI &mat_ids() const
    {
        return mat_IDs_;
    }

    bool is_fuel() const
    {
        return is_fuel_;
    }

private:
    // Pin ID
    const unsigned int id_;
    unsigned int mesh_id_;
    // Immutable reference to the pin mesh object (owned by CoreMesh)
    PinMesh const *pin_mesh_;
    // Material IDs to apply to each XS region of the pin mesh
    VecI mat_IDs_;
    //
    bool is_fuel_;
};

typedef std::shared_ptr<Pin> SP_Pin_t;
typedef std::unique_ptr<Pin> UP_Pin_t;

/**
 * \brief Given an XML node containing one or more \<pin\> tags, parse the
 * \ref Pin entries into a map.
 */
std::map<int, UP_Pin_t> ParsePins(const pugi::xml_node &input,
                                  const PinMesh_Map_t &meshes,
                                  const MaterialLib &mat_lib);
}
