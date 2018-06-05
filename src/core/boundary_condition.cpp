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

#include "core/boundary_condition.hpp"

#include "util/error.hpp"

namespace mocc {

BoundaryCondition::BoundaryCondition(int n_group,
                                     const AngularQuadrature &angquad,
                                     BC_Type_t bc, BC_Size_t n_bc,
                                     int dim)
    : BoundaryCondition(n_group, angquad, bc,
                        std::vector<BC_Size_t>(angquad.ndir(), n_bc), dim)
{
    return;
}

BoundaryCondition::BoundaryCondition(int n_group,
                                     const AngularQuadrature &angquad,
                                     BC_Type_t bc, std::vector<BC_Size_t> n_bc,
                                     int dim)
    : n_group_(n_group),
      n_angle_(n_bc.size()),
      bc_(bc),
      size_(n_bc),
      ang_quad_(angquad)
{
    assert((angquad.ndir() == (int)n_bc.size()) ||
           (angquad.ndir() / 2 == (int)n_bc.size()));

    assert((dim == 2) || (dim == 3));
    if (dim==2) {
        factor_ = 2.0;
    }
    else {
        factor_ = 1.0;
    }

    int n_angle = n_bc.size();

    bc_per_group_ = 0;
    for (auto n : n_bc) {
        bc_per_group_ += n[0] + n[1] + n[2];
    }

    size_t total_size = bc_per_group_ * n_group;
    data_.resize(total_size);
    offset_.resize(n_angle, 3);
    assert(offset_(0, blitz::Range::all()).isStorageContiguous());

    int offset = 0;
    int iang   = 0;
    for (auto n_ang : n_bc) {
        int iface = 0;
        for (auto n : n_ang) {
            offset_(iang, iface) = offset;
            offset += n;
            iface++;
        }
        iang++;
    }
    return;
}

BoundaryCondition::BoundaryCondition(const BoundaryCondition &rhs)
    : BoundaryCondition(rhs.n_group_, rhs.ang_quad_, rhs.bc_, rhs.size_,
            (rhs.factor_==2.0?2:3))
{
    return;
}

void BoundaryCondition::initialize_scalar(real_t val)
{
    // Start with all zeros
    data_ = 0.0;
    for (int group = 0; group < n_group_; group++) {
        // This is not a range-based loop, because n_angle_ is different
        // depending on the type of sweeper that made this
        // BoundaryCondition.
        for (int ang = 0; ang < n_angle_; ang++) {
            const auto &angle = ang_quad_[ang];
            for (auto norm : AllNormals) {
                Surface surf   = angle.upwind_surface(norm);
                auto face_pair = this->get_face(group, ang, norm);
                real_t *face   = face_pair.second;

                switch (bc_[(int)surf]) {
                case Boundary::VACUUM:
                    // Leave surface as all zeros
                    break;
                case Boundary::PRESCRIBED:
                    // Initialized to "1.0" for now, but will change it to
                    // be initialized in the Constructor via .h5 file. Then
                    // we will have to make sure this routine does not
                    // overwrite data_ to ZERO. It might make sense to have
                    // a separate initialization routine that does the
                    // initialization based on bc_.
                    for (int i = 0; i < face_pair.first; i++) {
                        face[i] = 1.0 / FPI * factor_;
                    }
                    break;
                case Boundary::PARALLEL:
                case Boundary::REFLECT:
                case Boundary::PERIODIC:
                    // initialize with the prescribed scalar
                    for (int i = 0; i < face_pair.first; i++) {
                        face[i] = val;
                    }
                    break;
                case Boundary::INVALID:
                    break;
                }
            }
        }
    }
}

void BoundaryCondition::initialize_spectrum(const ArrayB1 &spectrum)
{
    assert((int)spectrum.size() == n_group_);
    int it = 0;
    for (int ig = 0; ig < n_group_; ig++) {
        real_t val = spectrum(ig);
        data_(blitz::Range(it, it + bc_per_group_ - 1)) = val;
        it += bc_per_group_;
    }
    return;
}

void BoundaryCondition::update(int group, const BoundaryCondition &out)
{
    assert(out.n_group_ == 1);

    for (int iang = 0; iang < n_angle_; iang++) {
        this->update(group, iang, out);
    }
    return;
}

void BoundaryCondition::update(int group, int angle,
                               const BoundaryCondition &out)
{
    int group_offset = bc_per_group_ * group;

    for (Normal n : AllNormals) {
        int size    = size_[angle][(int)n];
        int iang_in = ang_quad_.reflect(angle, n);
        if (size == 0) {
            break;
        }
        assert(size == out.size_[iang_in][(int)n]);
        assert(iang_in < n_angle_);
        const auto &angle_in = ang_quad_[iang_in];
        int offset_in        = group_offset + offset_(iang_in, (int)n);
        int offset_out       = out.offset_(angle, (int)n);
        ;

        switch (bc_[(int)(angle_in.upwind_surface(n))]) {
        case Boundary::VACUUM:
            data_(blitz::Range(offset_in, offset_in + size - 1)) = 0.69780621255; //1.3956124251;
            //0.69780621255; //0.0;
            break;

        case Boundary::REFLECT:
            data_(blitz::Range(offset_in, offset_in + size - 1)) =
                out.data_(blitz::Range(offset_out, offset_out + size - 1));
            break;

        case Boundary::PRESCRIBED:
            break;

        default:
            throw EXCEPT("Unsupported boundary condition type");
        }
    }
    return;
}

std::ostream &operator<<(std::ostream &os, const BoundaryCondition &bc)
{
    os << "Boundary Condition:" << std::endl;
    for (auto b : bc.bc_) {
        os << b << std::endl;
    }
    os << std::endl;

    for (int igroup = 0; igroup < bc.n_group_; igroup++) {
        os << "Group: " << igroup << std::endl;

        for (int iang = 0; iang < bc.n_angle_; iang++) {
            os << "Angle: " << iang << std::endl;

            auto ang_size = bc.size_[iang];
            for (auto norm : AllNormals) {
                if (ang_size[(int)norm] > 0) {
                    os << norm << std::endl;
                    auto bvals = bc.get_face(igroup, iang, norm);
                    for (int i = 0; i < bvals.first; i++) {
                        os << bvals.second[i] << " ";
                    }
                    os << std::endl << std::endl;
                }
            }
        }
    }

    return os;
}
}
