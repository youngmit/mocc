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

#include "xs_mesh_region.hpp"

#include <algorithm>

namespace mocc {
    XSMeshRegion::XSMeshRegion( const VecI &fsrs,
            real_t *xstr,
            real_t *xsnf,
            real_t *xsch,
            real_t *xsf,
            real_t *xsrm,
            const ScatteringMatrix& xssc ):
        reg_(fsrs),
        xsmactr_(xstr),
        xsmacnf_(xsnf),
        xsmacf_(xsf),
        xsmacch_(xsch),
        xsmacrm_(xsrm),
        xsmacsc_(xssc)
    {
        is_fissile_ = false;
        for( int ig=0; ig<this->n_group(); ig++ ) {
            xsmacrm_[ig] = xsmactr_[ig] - xsmacsc_.self_scat(ig);
            if( xsmacnf_[ig] > 0.0 ) {
                is_fissile_ = true;
            }
        }
        return;
    }
} // namespace mocc
