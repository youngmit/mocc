#pragma once

#include "pugixml.hpp"
#include <blitz/array.h>

#include "mocc-core/angular_quadrature.hpp"
#include "mocc-core/global_config.hpp"

#include "sn/correction_data.hpp"
#include "sn/sn_sweeper_variant.hpp"
#include "sn/sn_sweeper_cdd.hpp"

#include "cmdo/moc_sweeper_2d3d.hpp"
#include "cmdo/source_2d3d.hpp"

namespace mocc {
    /**
     * This is an implementation of the 2D3D method. Each plane is treated with
     * a 2-D MoC sweeper, which produces the correction factors needed to treat
     * the entire system with a 3-D corrected diamond difference Sn sweeper.
     */
    class PlaneSweeper_2D3D: public TransportSweeper {
    public:
        PlaneSweeper_2D3D( const pugi::xml_node &input, const CoreMesh &mesh );

        void sweep( int group );

        void initialize();

        void get_pin_flux_1g( int ig, ArrayB1 &flux ) const;

        /**
         * Delegate to the subbordinate \ref sn::SnSweeper and \ref MoCSweeper.
         * Return the error from the MoC sweeper.
         */
        real_t set_pin_flux_1g( int group, const ArrayB1 &pin_flux ) {
            sn_sweeper_.set_pin_flux_1g( group, pin_flux );
            real_t diff = moc_sweeper_.set_pin_flux_1g( group, pin_flux );

            return diff;
        }

        void output( H5::CommonFG *file ) const;

        void homogenize( CoarseData &data ) const {
            throw EXCEPT("Not implemented");
        }

        /**
         * Associate the sweeper with a source. This has to do a little extra
         * work, since the Sn sweeper needs its own source.
         */
        virtual void assign_source( Source* source ) {
            assert( source != nullptr );
            source_ = source;
            moc_sweeper_.assign_source( source );
            /// \todo this static_cast is scary. Maybe think about relaxing the
            /// ownership of the source by FSS and allow the TS to figure out
            /// the types more explicitly...
            Source_2D3D *s = static_cast<Source_2D3D*>(source);
            sn_sweeper_.assign_source( s->get_sn_source() );
        }

        /**
         * Create a Source_2D3D object instead of the standard Source class.
         */
        UP_Source_t create_source() const {
            return UP_Source_t( new Source_2D3D( moc_sweeper_, sn_sweeper_ ) );
        }

        SP_XSMeshHomogenized_t get_homogenized_xsmesh() {
            return sn_sweeper_.get_homogenized_xsmesh();
        }

        /**
         * Override the default \ref TransportSweeper implementation to call the
         * method on one of the sub-sweepers.
         */
        void calc_fission_source( real_t k, ArrayF &fission_source ) const {
            moc_sweeper_.calc_fission_source( k, fission_source );
            return;
        }

        /**
         * \copybrief TransportSweeper::total_fission()
         * Override the default \ref TransportSweeper implementation to call the
         * method on one of the sub-sweepers.
         */
        real_t total_fission( bool old ) const {
            return sn_sweeper_.total_fission( old );
        }

        /**
         * Defer to the MoC and Sn sweepers.
         */
        void store_old_flux() {
            moc_sweeper_.store_old_flux();
            sn_sweeper_.store_old_flux();
            return;
        }

        void set_coarse_data( CoarseData *cd ) {
            coarse_data_ = cd;
            moc_sweeper_.set_coarse_data( cd );
            sn_sweeper_.set_coarse_data( cd );
        }

    private:
        // Parse the various options from the XML
        void parse_options( const pugi::xml_node &input );
        // Calculate transverse leakage based on the state of the coarse_data_
        // and apply to the MoC sweeper's source.
        void add_tl( int group );


        const CoreMesh& mesh_;
        sn::SnSweeperVariant<sn::CellWorker_CDD_DD> sn_sweeper_;
        MoCSweeper_2D3D moc_sweeper_;
        AngularQuadrature ang_quad_;
        std::shared_ptr<CorrectionData> corrections_;
        ArrayB2 tl_;

        // Sn-MoC residuals by group sweep
        std::vector<VecF> sn_resid_;

        // Options! Buttons and knobs!!!
        bool expose_sn_; // When asking the sweeper for pin flux, which one?
        bool do_snproject_;
        bool do_tl_;
        int n_inactive_moc_;
        int i_outer_;
    };
}
