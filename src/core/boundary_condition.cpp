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
