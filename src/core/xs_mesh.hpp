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

#include <vector>
#include "util/blitz_typedefs.hpp"
#include "util/fp_utils.hpp"
#include "util/global_config.hpp"
#include "core_mesh.hpp"
#include "output_interface.hpp"
#include "xs_mesh_region.hpp"

namespace mocc {
class XSMesh : public HasOutput {
public:
    // XSMesh provides its own facility to initialize itself from a \ref
    // CoreMesh
    XSMesh(const CoreMesh &mesh, MeshTreatment treatment);

    // Return the number of energy groups
    size_t n_group() const
    {
        return ng_;
    }

    // Iterators to the underlying vector
    const std::vector<XSMeshRegion>::const_iterator begin() const
    {
        return regions_.cbegin();
    }

    const std::vector<XSMeshRegion>::const_iterator end() const
    {
        return regions_.cend();
    }

    const XSMeshRegion &operator[](size_t i) const
    {
        assert(i >= 0 );
        assert(i < regions_.size());
        return regions_[i];
    }

    size_t size() const
    {
        return regions_.size();
    }

    const VecF &eubounds() const
    {
        return eubounds_;
    }

    /**
     * \brief Update macroscopic cross sections if needed
     *
     * Since the stock XSMesh only deals in un-homogenized, macroscopic
     * cross sections, this does nothing. When support for microscopic cross
     * sections is added, this will need to start doing some work.
     *
     * For right now, this is overridden in the \ref XSMeshHomogenized class
     * to calculate new homoginzed cross sections given a new state of the
     * FM scalar flux.
     */
    virtual void update()
    {
        // Do nothing for the regular XS Mesh... for now
        return;
    }

    virtual void output(H5Node &file) const
    {
        // Not really implementing for the general XS Mesh type.
        assert(false);
    }

    bool operator==(const XSMesh &other) const
    {
        if (regions_ != other.regions_) {
            return false;
        }
        return true;
    }

    bool operator!=(const XSMesh &other) const
    {
        return !(*this == other);
    }

    /**
     * \brief Return the number of regions that this XSMesh would expand into
     *
     * This is essentially the same as n_reg() for the sweeper with which the
     * XSMesh is associated.
     */
    int n_reg_expanded() const
    {
        return n_reg_expanded_;
    }

protected:
    /**
     * \brief Allocate space to store the actual cross sections.
     *
     * \param nxs the number of cross section mesh materials
     * \param ng the number of energy groups
     *
     * This is identical for all cross-section mesh types, so might as well
     * have it in one place.
     */
    void allocate_xs(int nxs, int ng)
    {
        auto shape = blitz::shape(nxs, ng);
        xstr_.resize(shape);
        xsnf_.resize(shape);
        xsch_.resize(shape);
        xsf_.resize(shape);
        xsrm_.resize(shape);
        auto test_slice(xstr_(0, blitz::Range::all()));
        assert(test_slice.isStorageContiguous());
    }

    size_t ng_;

    // Vector of xs mesh regions
    std::vector<XSMeshRegion> regions_;

    // Actual cross-section data
    ArrayB2 xstr_;
    ArrayB2 xsnf_;
    ArrayB2 xsch_;
    ArrayB2 xsf_;
    ArrayB2 xsrm_;

    // Energy group upper bounds
    VecF eubounds_;

    // Number of regions in the associated computational mesh
    int n_reg_expanded_;
};

typedef std::shared_ptr<XSMesh> SP_XSMesh_t;
}
