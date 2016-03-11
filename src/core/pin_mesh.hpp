#pragma once

#include <map>
#include <memory>

#include "pugixml.hpp"

#include "pin_mesh_base.hpp"
#include "pin_mesh_cyl.hpp"
#include "pin_mesh_rect.hpp"

/**
 * \page input Input Specification
 * \section pin_mesh <pin_mesh> Tag
 *
 * blah blah
 */

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
    PinMesh* PinMeshFactory(const pugi::xml_node &input);

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
     PinMesh_Map_t ParsePinMeshes( const pugi::xml_node &input );
}
