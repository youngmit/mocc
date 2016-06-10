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

#include "material_lib.hpp"

#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "error.hpp"
#include "files.hpp"
#include "material.hpp"

using std::stringstream;
using std::string;
using std::cout;
using std::cin;
using std::endl;

using std::regex;
using std::smatch;
using std::regex_match;

namespace mocc {

MaterialLib::MaterialLib()
{
    // do nothing
    return;
}

MaterialLib::MaterialLib(const pugi::xml_node &input) : n_material_(0)
{
    if (input.empty()) {
        throw EXCEPT("No material library specified.");
    }
    std::string matLibName = input.attribute("path").value();
    LogFile << "Using material library at: " << matLibName << std::endl;
    FileScrubber matLibFile;

    try {
        matLibFile = FileScrubber(matLibName.c_str(), "!");
    }
    catch (Exception e) {
        std::stringstream msg;
        std::cerr << e.what() << std::endl;
        msg << "Failed to open the cross-section library at: " << matLibName;
        throw EXCEPT(msg.str());
    }

    string line;
    // Read the first three lines to extract the library header
    // Description
    // Number of groups and number of materials
    {
        // skip the first line
        matLibFile.getline();
        stringstream inBuf(matLibFile.getline());
        inBuf >> n_grp_;
        if (inBuf.fail()) {
            throw EXCEPT("Failed to read number of groups!");
        }
        inBuf >> n_material_lib_;
        if (inBuf.fail()) {
            throw EXCEPT("Failed to read number of materials!");
        }
    }

    // Group boundaries
    {
        stringstream inBuf(matLibFile.getline());
        for (size_t i = 0; i < n_grp_; i++) {
            double bound;
            inBuf >> bound;
            g_bounds_.push_back(bound);
            if (inBuf.fail()) {
                throw EXCEPT("Trouble reading group bounds!");
            }
        }
    }

    // Read in material data
    for (size_t imat = 0; imat < n_material_lib_; imat++) {
        // Get the name of the material
        line = matLibFile.getline();

        regex headExp("^\\s*XSMACRO\\s+([^\\s]+)\\s+([0-9]+)\\s*$");
        smatch results;
        regex_match(line, results, headExp);
        string materialName = results[1].str();

        // Read in the non-scattering stuff
        VecF abs;
        VecF nuFiss;
        VecF fiss;
        VecF chi;
        for (size_t ig = 0; ig < n_grp_; ig++) {
            stringstream inBuf(matLibFile.getline());
            double val;

            inBuf >> val;
            abs.push_back(val);
            inBuf >> val;
            nuFiss.push_back(val);
            inBuf >> val;
            fiss.push_back(val);
            inBuf >> val;
            chi.push_back(val);
            if (inBuf.fail()) {
                throw EXCEPT("Trouble reading XS data from library!");
            }
        }

        // Read in the scattering table
        std::vector<VecF> scatTable;
        for (size_t ig = 0; ig < n_grp_; ig++) {
            stringstream inBuf(matLibFile.getline());
            VecF scatRow;
            for (size_t igg = 0; igg < n_grp_; igg++) {
                double val;
                inBuf >> val;
                scatRow.push_back(val);
            }
            scatTable.push_back(scatRow);
        }

        // produce a Material object and add it to the library
        lib_materials_.push_back(Material(abs, nuFiss, fiss, chi, scatTable));
        try {
            material_names_[materialName] = imat;
        }
        catch (...) {
            throw EXCEPT("Failed to add material from library. Duplicate "
                         "name?");
        }
    }

    // Parse material IDs
    for (auto mat = input.child("material"); mat;
         mat      = mat.next_sibling("material")) {
        this->assignID(mat.attribute("id").as_int(),
                       mat.attribute("name").value());
    }
    return;
}

MaterialLib::MaterialLib(FileScrubber &input) : n_material_(0)
{
    string line;
    // Read the first three lines to extract the library header
    // Description
    // Number of groups and number of materials
    {
        // skip the first line
        input.getline();
        stringstream inBuf(input.getline());
        inBuf >> n_grp_;
        if (inBuf.fail()) {
            throw EXCEPT("Failed to read number of groups!");
        }
        inBuf >> n_material_lib_;
        if (inBuf.fail()) {
            throw EXCEPT("Failed to read number of materials!");
        }
    }

    // Group boundaries
    {
        stringstream inBuf(input.getline());
        for (size_t i = 0; i < n_grp_; i++) {
            double bound;
            inBuf >> bound;
            g_bounds_.push_back(bound);
            if (inBuf.fail()) {
                throw EXCEPT("Trouble reading group bounds!");
            }
        }
    }

    // Read in material data
    for (size_t imat = 0; imat < n_material_lib_; imat++) {
        // Get the name of the material
        line = input.getline();

        regex headExp("^\\s*XSMACRO\\s+([^\\s]+)\\s+([0-9]+)\\s*$");
        smatch results;
        regex_match(line, results, headExp);
        string materialName = results[1].str();

        // Read in the non-scattering stuff
        VecF abs;
        VecF nuFiss;
        VecF fiss;
        VecF chi;
        for (size_t ig = 0; ig < n_grp_; ig++) {
            stringstream inBuf(input.getline());
            double val;

            inBuf >> val;
            abs.push_back(val);
            inBuf >> val;
            nuFiss.push_back(val);
            inBuf >> val;
            fiss.push_back(val);
            inBuf >> val;
            chi.push_back(val);
            if (inBuf.fail()) {
                throw EXCEPT("Trouble reading XS data from library!");
            }
        }

        // Read in the scattering table
        std::vector<VecF> scatTable;
        for (size_t ig = 0; ig < n_grp_; ig++) {
            stringstream inBuf(input.getline());
            VecF scatRow;
            for (size_t igg = 0; igg < n_grp_; igg++) {
                double val;
                inBuf >> val;
                scatRow.push_back(val);
            }
            scatTable.push_back(scatRow);
        }

        // produce a Material object and add it to the library
        lib_materials_.push_back(Material(abs, nuFiss, fiss, chi, scatTable));
        try {
            material_names_[materialName] = imat;
        }
        catch (...) {
            throw EXCEPT("Failed to add material from library. Duplicate "
                         "name?");
        }
    }
}

void MaterialLib::assignID(int id, std::string name)
{
    try {
        LogFile << "Mapping material '" << name << "' to ID " << id
                << std::endl;
        int mat_index = material_names_.at(name);
        assigned_materials_.push_back(lib_materials_[mat_index]);
        material_dense_index_[id] = n_material_;
        material_ids_[id]         = mat_index;
        n_material_++;
    }
    catch (std::out_of_range) {
        throw EXCEPT("Failed to map material to ID. Are you sure you "
                     "spelled it right?");
    }
    return;
}
}
