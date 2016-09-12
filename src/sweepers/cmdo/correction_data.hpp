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

#include <cassert>
#include <memory>
#include "util/blitz_typedefs.hpp"
#include "util/global_config.hpp"
#include "util/pugifwd.hpp"
#include "core/constants.hpp"
#include "core/core_mesh.hpp"
#include "core/output_interface.hpp"

namespace mocc {
/**
 * This class provides a storage scheme for the correction factors needed to
 * perform corrected diamond difference. The CDD Sn and MoC sweepers must be
 * provided with a reference to an object of this class to access and store
 * correction factors, respectively. Due to the relatively high
 * dimensionality of the data (space, angle, energy and cardinal direction
 * [X|Y]), instead of using a multidimensional array, we will instead use
 * accessor functions to get the data out of a dense linear representation.
 */
class CorrectionData : public HasOutput {
public:
    CorrectionData() : mesh_(nullptr), nreg_(0), nang_(0), ngroup_(0)
    {
        return;
    }

    CorrectionData(const CoreMesh &mesh, size_t nang, size_t ngroup)
        : mesh_(&mesh),
          nx_(mesh.nx()),
          ny_(mesh.ny()),
          nz_(mesh.macroplanes().size()),
          nreg_(nx_ * ny_ * nz_),
          nang_(nang),
          ngroup_(ngroup),
          alpha_(ngroup_, nang_, nreg_, 2),
          beta_(ngroup_, nang_, nreg_)
    {
        assert(alpha_.size() > 0);
        assert(beta_.size() > 0);
        assert(nx_ * ny_ * nz_ == nreg_);
        auto slice = alpha_(0, 0, blitz::Range::all(), blitz::Range::all());
        assert(slice.isStorageContiguous());

        alpha_ = 0.5;
        beta_  = 1.0;

        return;
    }

    ~CorrectionData()
    {
    }

    size_t size() const
    {
        return alpha_.size();
    }

    int n_cell() const
    {
        return nreg_;
    }

    inline real_t &alpha(int reg, int ang, int group, Normal norm)
    {
        return alpha_(group, ang, reg, (int)norm);
    }

    inline const real_t alpha(int reg, int ang, int group, Normal norm) const
    {
        return alpha_(group, ang, reg, (int)norm);
    }

    inline real_t &beta(int reg, int ang, int group)
    {
        return beta_(group, ang, reg);
    }

    inline const real_t beta(int reg, int ang, int group) const
    {
        return beta_(group, ang, reg);
    }

    /**
     * \brief Read correction factors from one or more HDF5 files, as
     * specified by \<data /\> tags
     */
    void from_data(const pugi::xml_node &input);

    void output(H5Node &file) const;

private:
    // Private methods to facilitate reading data from HDF5 files
    /**
     * \brief Read a single data file
     *
     * \param data an XML node containing a \<data\> specification
     *
     * This method reads data from a file specified in a \<data\> tag, applying
     * it to the entire problem geometry. This delegates to the \ref
     * read_data_single(const pugi::xml_node&, int, int) version, with the
     * macroplane bounds set to the extents of the core.
     */
    void read_data_single(const pugi::xml_node &data)
    {
        this->read_data_single(data, 0, mesh_->macroplanes().size() - 1);
        return;
    }

    /**
     * \brief Read a single data file
     *
     * \param data an XML node containing a \<data\> specification
     * \param bottom_plane the bottom macroplane index to apply the data to
     * \param top_plane the top macroplane index to apply the data to
     *
     * This method reads data from a file specified in a \<data\> tag, applying
     * it to the range of macroplanes specified.
     */
    void read_data_single(const pugi::xml_node &data, int bottom_plane,
                          int top_plane);

    /**
     * \brief Read data from all \<data\> tags in the passed XML node
     *
     * \param input an XML node that contains one or more \<data\> tag children
     *
     * This checks the \<data\> tags for validity, then delegates to \ref
     * read_data_single(const pugi::xml_node&, int, int)
     */
    void read_data_multi(const pugi::xml_node &input);
    const CoreMesh *mesh_;
    int nx_;
    int ny_;
    int nz_;
    int nreg_;
    int nang_;
    int ngroup_;

    ArrayB4 alpha_;
    ArrayB3 beta_;
};

typedef std::unique_ptr<CorrectionData> UP_CorrectionData_t;
}
