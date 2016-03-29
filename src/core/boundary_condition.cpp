#include "core/boundary_condition.hpp"

namespace mocc {
    std::ostream &operator<<(std::ostream &os, const BoundaryCondition &bc) {
        os << "Boundary Condition:" << std::endl;
        for(auto b: bc.bc_) {
            os << b << std::endl;
        }
        os << std::endl;

        for( int igroup=0; igroup<bc.n_group_; igroup++ ) {
            os << "Group: " << igroup << std::endl;

            for( int iang=0; iang<bc.n_angle_; iang++ ) {
                os << "Angle: " << iang << std::endl;

                auto ang_size = bc.size_[iang];
                for( auto norm: AllNormals ) {
                    if(ang_size[(int)norm] > 0) {
                        os << norm << std::endl;
                        auto bvals = bc.get_face(igroup, iang, norm);
                        for( int i=0; i<bvals.first; i++ ) {
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
