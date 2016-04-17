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

#include "xs_mesh_region.hpp"
#include "blitz_typedefs.hpp"
#include "core_mesh.hpp"
#include "fp_utils.hpp"
#include "global_config.hpp"
#include "output_interface.hpp"

namespace mocc {
    class XSMesh: public HasOutput {
    public:
        // Default constructor does almost nothing, and lets some other code
        // tell it what to do
        XSMesh() {}

        // XSMesh provides its own facility to initialize itself from a \ref
        // CoreMesh
        XSMesh( const CoreMesh& mesh );

        // Return the number of energy groups
        size_t n_group() const {
            return ng_;
        }

        // Iterators to the underlying vector
        const std::vector<XSMeshRegion>::const_iterator begin() const {
            return regions_.cbegin();
        }

        const std::vector<XSMeshRegion>::const_iterator end() const {
            return regions_.cend();
        }

        const XSMeshRegion& operator[]( size_t i ) const {
            return regions_[i];
        }

        size_t size() const {
            return regions_.size();
        }

        const VecF& eubounds() const {
            return eubounds_;
        }

        virtual void output( H5Node &file ) const {
            // Not really implementing for the general XS Mesh type.
            assert(false);
        }

        bool operator==( const XSMesh &other ) const {
            if( regions_ != other.regions_ ) {
                return false;
            }
            return true;
        }

        bool operator!=( const XSMesh &other ) const {
            return !(*this == other);
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
        void allocate_xs( int nxs, int ng ) {
            auto shape = blitz::shape(nxs, ng);
            xstr_.resize(shape);
            xsnf_.resize(shape);
            xsch_.resize(shape);
            xskf_.resize(shape);
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
        ArrayB2 xskf_;
        ArrayB2 xsrm_;

        // Energy group upper bounds
        VecF eubounds_;
    };

    typedef std::shared_ptr<XSMesh> SP_XSMesh_t;
}
