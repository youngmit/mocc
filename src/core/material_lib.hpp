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
#include <string>

#include "core/file_scrubber.hpp"
#include "core/global_config.hpp"
#include "core/material.hpp"
#include "core/pugifwd.hpp"

namespace mocc {
    typedef std::map<unsigned int, const Material*> MaterialMap;
    typedef std::vector<Material> MaterialVec;

    /**
     * The \ref MaterialLib stores a mapping of \ref Material objects, to be
     * used in constructing an \ref XSMesh.
     */
    class MaterialLib{
    public:
        MaterialLib();

        /**
         * Construct a \ref MaterialLib using a \ref FileScrubber configured to
         * parse an "MPACT user cross-section library." This whole approach is
         * admittedly wonky, but take a look at \ref CoreMesh::CoreMesh() to see
         * where this gets did.
         */
        MaterialLib(FileScrubber &input);

        MaterialLib( const pugi::xml_node &input );

        /**
         * Assign an ID to a material in the library.
         */
        void assignID(int id, std::string name);

        /**
         * Return the number of materials in the library
         */
        int n_materials() const {
            return n_material_;
        }

        /**
         * Return the map of materials by ID
         */
        const MaterialVec& materials() const {
            return lib_materials_;
        }

        /**
         * Return the index of the material given a material ID
         */
        int get_index_by_id( unsigned int id ) const {
            return material_dense_index_.at(id);
        }

        /**
         * Return a const reference to a material by ID
         */
        const Material& get_material_by_id( unsigned int id ) const {
            return lib_materials_[material_ids_.at(id)];
        }

        /**
         * Return a const refernce to the material indexed by ID
        */
        const Material& operator[]( int id ) const {
            return this->get_material_by_id( id );
        }

        /**
         * Return the number of groups spanned by the library
         */
        int n_group() const {
            return n_grp_;
        }

        /**
         * Return the group bounds
         */
        const VecF& g_bounds() const {
            return g_bounds_;
        }

        MaterialVec::const_iterator begin() const {
            return assigned_materials_.cbegin();
        }

        MaterialVec::const_iterator end() const {
            return assigned_materials_.cend();
        }

        /**
         * \brief Return whether the passed material ID is defined
         */
        bool has( int id ) const {
            return material_ids_.count(id) > 0;
        }

    private:
        // Vector storing all of the materials in the library.
        MaterialVec lib_materials_;

        // Vector of actual assigned materials. A shame to do a copy, but avoids
        // having to do all sorts of iterator jiggery later
        MaterialVec assigned_materials_;

        // Map from a material name to its corresponding index in materials_
        std::map<std::string, int> material_names_;

        // Map from a material ID to the corresponding index in the vector of
        // materials
        std::map<int, int> material_ids_;

        // Map from a material ID to a corresponding index in a dense index
        // space ( 0, n_material() )
        std::map<int, int> material_dense_index_;

        // Number of energy groups for which all materials in the library are
        // defined
        unsigned int n_grp_;

        // Number of materials that have been mapped to an ID
        unsigned int n_material_;

        // Number of materials present in the library itself (>= n_material_)
        unsigned int n_material_lib_;

        // Upper bound for each of the energy groups
        VecF g_bounds_;

        // Descriptive string for the material library
        std::string m_description;
    };
}
