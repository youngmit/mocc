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

#include "scattering_matrix.hpp"

#include <iomanip>
#include <iostream>

#include "error.hpp"

namespace {
using namespace mocc;
ArrayB2 convert_vvec(const std::vector<VecF> &scat)
{
    // Convert the vector of vectors to and Array2, then call that version
    // First make sure its sqare
    for (const auto &row : scat) {
        if (row.size() != scat.size()) {
            throw EXCEPT("Scattering matrix input is not square");
        }
    }

    ArrayB2 array(scat.size(), scat.size(), 0.0);
    int irow = 0;
    for (const auto &row : scat) {
        int icol = 0;
        for (const auto &v : row) {
            array(irow, icol) = v;
            icol++;
        }
        irow++;
    }

    return array;
}
}

namespace mocc {
ScatteringMatrix::ScatteringMatrix(const std::vector<VecF> &scat)
    : ScatteringMatrix(convert_vvec(scat))
{
    return;
}

ScatteringMatrix::ScatteringMatrix(const ArrayB2 &scat)
{
    // Make sure that scat is square
    assert(scat.extent(0) == scat.extent(1));

    // Imply ng_ from the size of the passed-in vectors
    ng_  = scat.extent(0);
    out_ = VecF(ng_, 0.0);
    rows_.reserve(ng_);

    // Determien the size of scat_ and save the bounds for ScatteringRows
    int size = 0;
    std::vector<std::pair<int, int>> bounds;
    bounds.reserve(ng_);

    for (int to = 0; to < ng_; to++) {
        std::pair<int, int> these_bounds;
        for (int from = 0; from < ng_; from++) {
            these_bounds.first = from;
            if (scat(to, from) > 0.0) {
                break;
            }
        }
        for (int from = ng_ - 1; from >= 0; from--) {
            these_bounds.second = from;
            if (scat(to, from) > 0.0) {
                break;
            }
        }
        // Handle empty scattering rows. If the min and max g end up on the
        // other side of the group range, we didn't find any non-zero cross
        // sections. Set min_g and max_g to the current group. When we pack
        // the cross sections into their dense representations, we make sure
        // to add a corresponding zero
        if (these_bounds.first == ng_ - 1 && these_bounds.second == 0) {
            these_bounds.first  = to;
            these_bounds.second = to;
        }

        // Also make sure that the scattering rows include the self-scatter
        // entry
        if (these_bounds.first > to) {
            these_bounds.first = to;
        }
        if (these_bounds.second < to) {
            these_bounds.second = to;
        }
        size += these_bounds.second - these_bounds.first + 1;
        bounds.push_back(these_bounds);
    }

    scat_.reserve(size);

    int to  = 0;
    int pos = 0;
    for (auto &these_bounds : bounds) {
        for (int from = these_bounds.first; from <= these_bounds.second;
             from++) {
            scat_.push_back(scat(to, from));
            out_[from] += scat(to, from);
        }
        rows_.push_back(ScatteringRow(these_bounds.first, these_bounds.second,
                                      &scat_[pos]));
        to++;
        pos += these_bounds.second - these_bounds.first + 1;
    }
}

std::ostream &operator<<(std::ostream &os, const ScatteringMatrix &scat_mat)
{
    for (auto &row : scat_mat.rows_) {
        int gmin = row.min_g;
        int gmax = row.max_g;
        for (int ig = 0; ig < gmin; ig++) {
            os << std::setw(12) << 0.0;
        }

        for (int ig = gmin; ig <= gmax; ig++) {
            os << std::setw(12) << row[ig];
        }

        for (int ig = gmax + 1; ig < scat_mat.ng_; ig++) {
            os << std::setw(12) << 0.0;
        }
        os << std::endl;
    }
    return os;
}
}
