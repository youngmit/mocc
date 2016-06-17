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

#include "core/blitz_typedefs.hpp"
#include "core/constants.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"
#include "core/output_interface.hpp"
#include "core/pugifwd.hpp"

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
        CorrectionData( ):
            mesh_( nullptr ),
            nreg_( 0 ),
            nang_( 0 ),
            ngroup_( 0 )
        {
            return;
        }

        CorrectionData( const Mesh &mesh, size_t nang, size_t ngroup ):
            mesh_( &mesh ),
            nreg_( mesh.n_pin() ),
            nx_( mesh.nx() ),
            ny_( mesh.ny() ),
            nz_( mesh.nz() ),
            nang_( nang ),
            ngroup_( ngroup ),
            alpha_( ngroup_, nang_, nreg_, 2 ),
            beta_( ngroup_, nang_, nreg_ )
        {
            assert( alpha_.size() > 0 );
            assert( beta_.size() > 0 );
            assert( nx_*ny_*nz_ == nreg_ );
            auto slice = alpha_(0, 0, blitz::Range::all(), blitz::Range::all());
            assert( slice.isStorageContiguous() );

            alpha_ = 0.5;
            beta_ = 1.0;

            return;
        }

        ~CorrectionData() { }

        size_t size() {
            return alpha_.size();
        }


        inline real_t& alpha( int reg, int ang, int group, Normal norm )
        {
            return alpha_( group, ang, reg, (int)norm );
        }

        inline const real_t alpha( int reg, int ang, int group,
                Normal norm ) const
        {
            return alpha_( group, ang, reg, (int)norm );
        }

        inline real_t& beta( int reg, int ang, int group ) {
            return beta_( group, ang, reg );
        }

        inline const real_t beta( int reg, int ang, int group ) const
        {
            return beta_( group, ang, reg );
        }

        /**
         * \brief Read correction factors from one or more HDF5 files, as
         * specified by \<data /\> tags
         */
        void from_data( const pugi::xml_node &input );

        void output( H5Node &file ) const;

    private:
        const Mesh *mesh_;
        int nreg_;
        int nx_;
        int ny_;
        int nz_;
        int nang_;
        int ngroup_;

        ArrayB4 alpha_;
        ArrayB3 beta_;
    };

    typedef std::unique_ptr<CorrectionData> UP_CorrectionData_t;
}
