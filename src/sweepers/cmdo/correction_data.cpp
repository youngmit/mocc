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
    // First, validate the data tags. Make sure that they are the right
    // size and have cover all planes in the mesh.
    {
        std::vector<bool> plane_data(nz_, false);
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

            if ((bot_plane < 0) || (bot_plane >= (int)mesh_->nz())) {
                throw EXCEPT("Invalid bottom_plane");
            }
            if ((top_plane < 0) || (top_plane >= (int)mesh_->nz())) {
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
                   "planes:"
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
        H5Node h5f(data.attribute("file").value(), H5Access::READ);

        for (int ig = 0; ig < ngroup_; ig++) {
            for (int iang = 0; iang < nang_; iang++) {
                // Gobble that data. Om nom nom!
                {
                    std::stringstream path;
                    path << "/alpha_x/" << std::setfill('0') << std::setw(3)
                         << ig << "_" << std::setfill('0') << std::setw(3)
                         << iang;

                    h5f.read_1d(path.str(), slice);
                    assert(slice.size() == mesh_->n_cell_plane());
                    for (int ip = bot_plane; ip <= top_plane; ip++) {
                        int stt = mesh_->plane_cell_begin(ip);
                        int stp = mesh_->plane_cell_end(ip) - 1;
                        alpha_(ig, iang, blitz::Range(stt, stp),
                               (int)Normal::X_NORM) = slice;
                    }
                }

                {
                    std::stringstream path;
                    path << "/alpha_y/" << std::setfill('0') << std::setw(3)
                         << ig << "_" << std::setfill('0') << std::setw(3)
                         << iang;
                    h5f.read_1d(path.str(), slice);
                    for (int ip = bot_plane; ip <= top_plane; ip++) {
                        int stt = mesh_->plane_cell_begin(ip);
                        int stp = mesh_->plane_cell_end(ip) - 1;
                        alpha_(ig, iang, blitz::Range(stt, stp),
                               (int)Normal::Y_NORM) = slice;
                    }
                }
                {
                    std::stringstream path;
                    path << "/beta/" << std::setfill('0') << std::setw(3) << ig
                         << "_" << std::setfill('0') << std::setw(3) << iang;
                    h5f.read_1d(path.str(), slice);
                    for (int ip = bot_plane; ip <= top_plane; ip++) {
                        int stt = mesh_->plane_cell_begin(ip);
                        int stp = mesh_->plane_cell_end(ip) - 1;
                        beta_(ig, iang, blitz::Range(stt, stp)) = slice;
                    }
                }
            }
        }
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
