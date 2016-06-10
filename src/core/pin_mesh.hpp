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

#include "pin_mesh_base.hpp"
#include "pin_mesh_cyl.hpp"
#include "pin_mesh_rect.hpp"
#include "pugifwd.hpp"

namespace mocc {
typedef std::shared_ptr<PinMesh> SP_PinMesh_t;
typedef std::unique_ptr<PinMesh> UP_PinMesh_t;
typedef std::map<int, UP_PinMesh_t> PinMesh_Map_t;

/**
 * \brief Construct a pin mesh, and return a \c unique_ptr to it.
 *
 * \param input an XML node, which should be a valid \c \<mesh\> tag
 *
 * Given XML input, determine the concrete type of \ref PinMesh to make,
 * construct it, and pass a managed pointer to it back to the caller.
 */
PinMesh *PinMeshFactory(const pugi::xml_node &input);

/**
 * \brief Parse all \ref PinMesh es in a given XML node, and return a map of
 * pointers, keyed by their respective IDs.
 *
 * \param input a reference to an XML node containing one or more \<mesh\>
 * children.
 *
 * This walks through the passed XML node, parsing each pin mesh specified
 * within.
 */
PinMesh_Map_t ParsePinMeshes(const pugi::xml_node &input);
}
