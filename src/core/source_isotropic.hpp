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

#include "core/source.hpp"

namespace mocc {
    /**
     * This extends the virtual \ref Source type for use as an isotropic source
     * for MoC sweepers.
     */
    class SourceIsotropic : public Source {
    public:
        SourceIsotropic( int nreg, const XSMesh* xs_mesh, const ArrayB2& flux ):
            Source( nreg, xs_mesh, flux),
            q_( nreg )
        {
            q_.fill(0.0);
            return;
        }

        virtual void self_scatter( size_t ig,
                const ArrayB1 &xstr = ArrayB1(0) );

        const VectorX& get_transport( int iang ) const {
            return q_;
        }

    protected:
        // The source, including self-scatter. This is stored separately from
        // source_1g_ so that the self_scatter method may be called multiple
        // times without having to completely reconstruct the source. All calls
        // to get_transport() will return a reference to this vector.
        VectorX q_;
    };
}
