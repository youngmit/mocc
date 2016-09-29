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

#include "correction_data.hpp"

#include <iomanip>
#include "pugixml.hpp"
#include "util/files.hpp"
#include "util/string_utils.hpp"

namespace mocc {
void CorrectionData::from_data(const pugi::xml_node &input)
{
    if (input.child("data").empty()) {
        // There isn't actually any data. Go ahead and return
        LogFile << "CorrectionData::from_data() was called, but <data/> "
                   "was empty or non-existent.";
        return;
    }

    LogFile << "Loading CDD data from file(s)." << std::endl;

    int np = mesh_->macroplanes().size();

    int n_data = 0;

    for (const auto &d : input.children("data")) {
        n_data++;
        LogFile << "Looking for correction data in file: "
                << d.attribute("file").value() << std::endl;
    }

    bool single_file = false;
    if (n_data == 1) {
        // peek at the size of the data and see if it matches the entire geom
        H5Node h5d(input.child("data").attribute("file").value(),
                   H5Access::READ);
        auto dims = h5d.dimensions("alpha_x/000/000");
        if ((dims[2] == mesh_->nx()) && (dims[1] == mesh_->ny()) &&
            ((int)dims[0] == np)) {
            single_file = true;
        }
    }

    if (single_file) {
        this->read_data_single(input.child("data"));
    } else {
        this->read_data_multi(input);
    }

    return;
}

void CorrectionData::read_data_single(const pugi::xml_node &data,
                                      int bottom_plane, int top_plane)
{
    // Make sure that if specified, the [macro]plane bounds correspond to the
    // mesh
    {
        int bot_plane_in = 0;
        int top_plane_in = mesh_->macroplanes().size() - 1;
        if (!data.attribute("top_plane").empty()) {
            top_plane_in = data.attribute("top_plane").as_int(-1);
        }
        if (!data.attribute("bottom_plane").empty()) {
            bot_plane_in = data.attribute("bottom_plane").as_int(-1);
        }

        if ((bot_plane_in != bottom_plane) || (top_plane_in != top_plane)) {
            throw EXCEPT("Plane bounds specified, but invalid");
        }
    }

    H5Node h5d(data.attribute("file").value(), H5Access::READ);

    // allocate a single buffer to read our data into for each group/angle
    ArrayB1 inbuf(nx_ * ny_);

    for (int ig = 0; ig < ngroup_; ig++) {
        for (int iang = 0; iang < nang_; iang++) {
            // alpha x
            {
                std::stringstream path;
                path << "/alpha_x/" << std::setw(3) << std::setfill('0') << ig
                     << "/" << std::setw(3) << std::setfill('0') << iang;
                try {
                    h5d.read(path.str(), inbuf);
                } catch (Exception e) {
                    throw EXCEPT_E("Failed to read alpha x", e);
                }
                for (int ip = bottom_plane; ip <= top_plane; ip++) {
                    int stt = mesh_->plane_cell_begin(ip);
                    int stp = mesh_->plane_cell_end(ip) - 1;
                    alpha_(ig, iang, blitz::Range(stt, stp),
                           (int)Normal::X_NORM) = inbuf;
                }
            }
            // alpha y
            {
                std::stringstream path;
                path << "/alpha_y/" << std::setw(3) << std::setfill('0') << ig
                     << "/" << std::setw(3) << std::setfill('0') << iang;
                try {
                    h5d.read(path.str(), inbuf);
                } catch (Exception e) {
                    throw EXCEPT_E("Failed to read alpha y", e);
                }
                for (int ip = bottom_plane; ip <= top_plane; ip++) {
                    int stt = mesh_->plane_cell_begin(ip);
                    int stp = mesh_->plane_cell_end(ip) - 1;
                    alpha_(ig, iang, blitz::Range(stt, stp),
                           (int)Normal::Y_NORM) = inbuf;
                }
            }
            // beta
            {
                std::stringstream path;
                path << "/beta/" << std::setw(3) << std::setfill('0') << ig
                     << "/" << std::setw(3) << std::setfill('0') << iang;
                try {
                    h5d.read(path.str(), inbuf);
                } catch (Exception e) {
                    throw EXCEPT_E("Failed to read beta", e);
                }
                for (int ip = bottom_plane; ip <= top_plane; ip++) {
                    int stt = mesh_->plane_cell_begin(ip);
                    int stp = mesh_->plane_cell_end(ip) - 1;
                    beta_(ig, iang, blitz::Range(stt, stp)) = inbuf;
                }
            }
        }
    }

    return;
}

void CorrectionData::read_data_multi(const pugi::xml_node &input)
{
    // First, validate the data tags. Make sure that they are the right
    // size and cover all planes in the mesh.
    int np = mesh_->macroplanes().size();
    {
        std::vector<bool> plane_data(np, false);
        for (auto data = input.child("data"); data;
             data      = data.next_sibling("data")) {
            int top_plane = np - 1;
            int bot_plane = 0;
            if (!data.attribute("top_plane").empty()) {
                top_plane = data.attribute("top_plane").as_int();
            }
            if (!data.attribute("bottom_plane").empty()) {
                bot_plane = data.attribute("bottom_plane").as_int();
            }

            if ((bot_plane < 0) || (bot_plane >= np)) {
                throw EXCEPT("Invalid bottom_plane");
            }
            if ((top_plane < 0) || (top_plane >= np)) {
                throw EXCEPT("Invalid top_plane");
            }

            if (data.attribute("file").empty()) {
                throw EXCEPT("No file specified.");
            }

            for (int ip = bot_plane; ip <= top_plane; ip++) {
                if (plane_data[ip]) {
                    std::stringstream msg;
                    msg << "Plane data is over-specified. Look at plane " << ip;
                    throw EXCEPT(msg.str());
                }
                plane_data[ip] = true;
            }
        }
        LogFile << "Correction data is being specified for the following "
                   "macroplanes:"
                << std::endl;
        LogFile << print_range(plane_data) << std::endl;
    }

    ArrayB1 slice(mesh_->n_cell_plane());

    for (auto data = input.child("data"); data;
         data      = data.next_sibling("data")) {
        int top_plane = 0;
        int bot_plane = 0;
        if (!data.attribute("top_plane").empty()) {
            top_plane = data.attribute("top_plane").as_int();
        }
        if (!data.attribute("bottom_plane").empty()) {
            bot_plane = data.attribute("bottom_plane").as_int();
        }

        this->read_data_single(data, bot_plane, top_plane);
    }
    return;
}

void CorrectionData::output(H5Node &file) const
{
    VecI dims;
    dims.push_back(nz_);
    dims.push_back(ny_);
    dims.push_back(nx_);
    int n = nx_ * ny_ * nz_;

    file.create_group("alpha_x");
    file.create_group("alpha_y");
    file.create_group("beta");

    // Declare slice storage
    ArrayB1 slice(n);

    for (int g = 0; g < ngroup_; g++) {
        std::stringstream path;
        path << "alpha_x/" << std::setfill('0') << std::setw(3) << g;
        auto ax_g = file.create_group(path.str());

        path.str("");
        path << "alpha_y/" << std::setfill('0') << std::setw(3) << g;
        auto ay_g = file.create_group(path.str());

        path.str("");
        path << "beta/" << std::setfill('0') << std::setw(3) << g;
        auto beta_g = file.create_group(path.str());

        for (int a = 0; a < nang_; a++) {
            {
                slice = beta_(g, a, blitz::Range::all());
                std::stringstream setname;
                setname << std::setfill('0') << std::setw(3) << a;
                beta_g.write(setname.str(), slice, dims);
            }

            {
                slice = alpha_(g, a, blitz::Range::all(), (int)Normal::X_NORM);
                std::stringstream setname;
                setname << std::setfill('0') << std::setw(3) << a;
                ax_g.write(setname.str(), slice, dims);
            }

            {
                slice = alpha_(g, a, blitz::Range::all(), (int)Normal::Y_NORM);
                std::stringstream setname;
                setname << std::setfill('0') << std::setw(3) << a;
                ay_g.write(setname.str(), slice, dims);
            }
        }
    }
    return;
}
}
