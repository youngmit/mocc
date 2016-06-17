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

#include "core/global_config.hpp"

#include "sweepers/moc/moc_sweeper.hpp"

#include "correction_data.hpp"

namespace mocc { namespace cmdo {
    class MoCSweeper_2D3D: public moc::MoCSweeper {
    public:
        MoCSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        /**
         * \brief Assign correction and cross-section coupling.
         */
        void set_coupling( std::shared_ptr<CorrectionData> data,
                SP_XSMeshHomogenized_t xsmesh ) {
            if( corrections_ || sn_xs_mesh_ ) {
                throw EXCEPT( "Correction data already assigned." );
            }
            corrections_ = data;
            sn_xs_mesh_ = xsmesh;
        }

        /**
         * \brief Allocate space internally to store coupling coefficients and
         * cross sections. Mainly useful for one-way coupling.
         */
        void set_self_coupling() {
            internal_coupling_ = true;
            corrections_ = std::shared_ptr<CorrectionData>(
                    new CorrectionData( mesh_,
                    ang_quad_.ndir()/2, xs_mesh_->n_group() ) );
            sn_xs_mesh_ = this->get_homogenized_xsmesh();
            sn_xs_mesh_->set_flux( flux_ );
        }

        /**
         * \brief \copybrief moc::MoCSweeper::expand_xstr()
         *
         * \copydetails moc::MoCSweeper::expand_xstr()
         *
         * In addition to what is performed by the \ref
         * moc::MoCSweeper::expand_xstr() method, this also expands the
         * unaltered transport cross sections to a separate array for the \ref
         * CurrentCorrections. We do not want to use altered cross section for
         * the Sn sweeper, which treats the axial dimension explicitly.
         */
        void expand_xstr( int group );

        /**
         * \brief Extend output() to export correction factors and homogenized
         * cross sections if the sweeper is internally coupled. Again, this is
         * only for the one-way coupling case to output data.
         */
        void output( H5Node &node ) const;

    private:
        void sweep1g_final( int group );

        /**
         * Given homogenized angular flux and total cross section data,
         * calculate the correction factors for CDD
         */
        void calculate_corrections( size_t ang, size_t group, ArrayF flux_surf,
                ArrayF flux_node, ArrayF sigt );

        std::shared_ptr<CorrectionData> corrections_;

        /**
         * The transport cross sections for the current group, unaltered due to
         * source splitting.
         */
        ArrayB1 xstr_true_;

        std::shared_ptr<XSMeshHomogenized> sn_xs_mesh_;

        bool internal_coupling_;

        std::vector<std::vector<std::array<real_t, 3>>> correction_residuals_;
    };
} } // Namespace mocc::cmdo
