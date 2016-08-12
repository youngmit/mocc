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
#include "util/global_config.hpp"
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

    /**
     * \brief Return the total surface area of the Pin
     *
     * \note Pins possess no concept of axial height, so their "volume" is
     * really surface area
     */
    real_t area() const
    {
        return pin_mesh_->area();
    }

    /**
    * \brief Return a vector of the areas of each mesh region
    *
    * \note Pins possess no concept of axial height, so their "volume" is
    * really surface area
    */
    const VecF &areas() const
    {
        return pin_mesh_->areas();
    }

    const VecI &mat_ids() const
    {
        return mat_IDs_;
    }

    bool is_fuel() const
    {
        return is_fuel_;
    }

    /**
     * \brief Return whether two pins are the same
     *
     * This is just checking the ids of the pins. At some point it may be useful
     * to check other porperties of the Pin since it is possible to want to
     * compare pins from different index universes, but for now this is
     * sufficient
     */
    bool operator==(const Pin &other) const
    {
        return id_ == other.id_;
    }

    friend std::ostream &operator<<(std::ostream &os, const Pin &pin);

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
