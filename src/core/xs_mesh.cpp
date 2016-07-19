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

#include "xs_mesh.hpp"

#include <iostream>
#include <map>
#include <memory>
#include "util/blitz_typedefs.hpp"
#include "util/files.hpp"
#include "util/global_config.hpp"

namespace mocc {
XSMesh::XSMesh(const CoreMesh &mesh, MeshTreatment treatment): state_(0)
{

    LogFile << "Initializing XS Mesh... ";

    const MaterialLib &mat_lib = mesh.mat_lib();

    // Assume the same number of groups as the source material library
    ng_ = mat_lib.n_group();

    // Get energy group bounds
    eubounds_ = mat_lib.g_bounds();

    // Set up the XS mesh regions. This essentially boils down to generating a
    // map from material index to the flat source region indices that are filled
    // by the indexed material. After that, everything else should be quite
    // similar.
    std::map<int, VecI> fsr_map;

    if (treatment == MeshTreatment::TRUE) {
        int ireg = 0;
        for (const auto &pini : mesh) {
            const PinMesh &pm   = pini->mesh();
            const VecI &mat_ids = pini->mat_ids();
            int ixsreg          = 0;
            for (auto &mat_id : mat_ids) {
                for (size_t reg = 0; reg < pm.n_fsrs(ixsreg); reg++) {
                    fsr_map[mat_id].push_back(ireg);
                    ireg++;
                }
                ixsreg++;
            }
        }
    } else if (treatment == MeshTreatment::PLANE) {
        int ireg = 0;
        int iz   = 0;
        for (const auto &block : mesh.subplane()) {
            const Plane &plane = mesh.unique_plane(mesh.unique_plane_id(iz));
            for (const auto &lattice : plane) {
                for (const auto &pin : *lattice) {
                    const VecI &mat_ids = pin->mat_ids();
                    const PinMesh &pm = pin->mesh();
                    int ixsreg        = 0;
                    for (const auto &mat_id : mat_ids) {
                        for(unsigned reg=0; reg<pm.n_fsrs(ixsreg); ++reg ){
                            assert(ireg < nreg);
                            fsr_map[mat_id].push_back(ireg);
                            ireg++;
                        }
                        ixsreg++;
                    }
                }
            }
            iz += block;
        }
    } else {
        // Should be using the homogenized class. This would be nice to merge at
        // some point though
        return;
    }

    int n_xsreg = fsr_map.size();

    this->allocate_xs(n_xsreg, ng_);

    // The ids/keys in fsr_map correspond to the user-specified IDs in the
    // material library. We want to cast this into a contiguous, zero-based
    // index space for internal storage and saner indexing, hence the imat
    // counter.
    n_reg_expanded_ = 0;
    VecI mat_ids(n_xsreg);
    int imat = 0;
    for (const auto &mat_pair : fsr_map) {
        mat_ids[imat]   = mat_pair.first;
        const auto &mat = mat_lib[mat_pair.first];
        xstr_(imat, blitz::Range::all()) = mat.xstr();
        xsnf_(imat, blitz::Range::all()) = mat.xsnf();
        xsch_(imat, blitz::Range::all()) = mat.xsch();
        xsf_(imat, blitz::Range::all())  = mat.xsf();
        n_reg_expanded_ += mat_pair.second.size();

        imat++;
    }

    // Preallocate space for the regions. Saves on lots of copies for large
    // xsmeshes.
    regions_.reserve(fsr_map.size());
    imat = 0;
    for (auto &reg_pair : fsr_map) {
        const auto &mat = mat_lib[mat_ids[imat]];
        regions_.emplace_back(reg_pair.second, &xstr_(imat, 0), &xsnf_(imat, 0),
                              &xsch_(imat, 0), &xsf_(imat, 0), &xsrm_(imat, 0),
                              mat.xssc());
        imat++;
    }

    LogFile << "done." << std::endl;

    return;
}

std::ostream &operator<<(std::ostream &os, const XSMeshRegion &xsr)
{
    int ng = xsr.xsmacsc_.n_group();
    os << "Transport: " << std::endl;
    for (int ig = 0; ig < ng; ig++) {
        os << xsr.xsmactr_[ig];
    }
    os << std::endl;

    os << "nu-fission: " << std::endl;
    for (int ig = 0; ig < ng; ig++) {
        os << xsr.xsmacnf_[ig];
    }
    os << std::endl;

    os << "chi: " << std::endl;
    for (int ig = 0; ig < ng; ig++) {
        os << xsr.xsmacch_[ig];
    }
    os << std::endl;

    os << "removal: " << std::endl;
    for (int ig = 0; ig < ng; ig++) {
        os << xsr.xsmacrm_[ig];
    }
    os << std::endl;

    os << "Scattering matrix:" << std::endl;
    os << xsr.xsmacsc_ << std::endl;

    return os;
}
}
