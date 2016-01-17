#pragma once

#include "core/global_config.hpp"

#include "sweepers/moc/moc_sweeper.hpp"

#include "sweepers/sn/correction_data.hpp"

namespace mocc {
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
         * \brief Extend output() to export correction factors and homogenized
         * cross sections if the sweeper is internally coupled. Again, this is
         * only for the one-way coupling case to output data.
         */
        void output( H5::CommonFG *node ) const;

    private:
        void sweep1g_final( int group );

        /**
         * Given homogenized angular flux and total cross section data,
         * calculate the correction factors for CDD
         */
        void calculate_corrections( size_t ang, size_t group, ArrayF flux_surf,
                ArrayF flux_node, ArrayF sigt );

        std::shared_ptr<CorrectionData> corrections_;

        std::shared_ptr<XSMeshHomogenized> sn_xs_mesh_;

        bool internal_coupling_;
    };
}
