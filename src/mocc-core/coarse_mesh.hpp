#pragma once

#include <string>
#include <map>

#include "global_config.hpp"
#include "mesh.hpp"
#include "core_mesh.hpp"

namespace mocc {
    /**
    * The CoarseMesh type is sort of a one-stop shop for the data that might be
    * necessary for solvers that operate on pin-homogenized meshes. It doesnt
    * really store that much in the way of mesh properties since we guarantee
    * that it conforms to the CoreMesh, and is therefore a simple structured
    * grid. 
    *
    * The CoarseMesh provides indexing for both the volumetric cells
    * (correspontant to the pin cells on the CoreMesj), as well as the interface
    * surfaces between the pin cells and the boundary of the problem domain.
    * These indices, while straightforward, should be obtained using the
    * surface_index() and cell_index() methods. See the documentaion for those
    * methods to see how the indexing is performed. 
    *
    * \todo At some point it might be nice to extend the base Mesh type and use
    * this as the Sn mesh, too.
    */
    class CoarseMesh {
    public:
        CoarseMesh( const CoreMesh& mesh ):
            mesh_(mesh)    
        { return; };

        unsigned int cell_index() {
            return 0;
        }

        unsigned int surface_index( int cell, Surface surf ) {
            return 0;
        }
    private:
        /**
         * Const reference to a CoreMesh. The CoarseMesh will be conformant to
         * that mesh
        */
        const CoreMesh& mesh_;

        /**
        * The actual data. Since what needs to be stored is problem-dependent,
        * we use a map, keyed on strings to keep things dynamic.
        */
        std::map<std::string, VecF> data_;
    };
}
