#include "plane_sweeper_2d3d.hpp"

namespace mocc {
    PlaneSweeper_2D3D::PlaneSweeper_2D3D( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        sn_sweeper_(input.child("sn_sweeper"), mesh),
        moc_sweeper_(input.child("moc_sweeper"), mesh),
        ang_quad_(moc_sweeper_.get_ang_quad())
    {
        /// \todo initialize the rest of the data members on the TS base type
        xs_mesh_ = moc_sweeper_.get_xs_mesh();
        n_reg_ = moc_sweeper_.n_reg();
        ng_ = xs_mesh_->n_group();

        corrections_ = CorrectionData( sn_sweeper_.n_reg(), ang_quad_.ndir()/2,
                ng_ );

        sn_sweeper_.set_corrections( &corrections_ );
        moc_sweeper_.set_corrections( &corrections_ );
        
        return;
    }

    void PlaneSweeper_2D3D::sweep( int group ) {
        moc_sweeper_.sweep( group );
        sn_sweeper_.sweep( group );
        
    }

    void PlaneSweeper_2D3D::initialize() {
        sn_sweeper_.initialize();
        moc_sweeper_.initialize();
    }

    void PlaneSweeper_2D3D::get_pin_flux( int ig, VecF &flux ) const {
        /// \todo implement. probably just call it on the Sn sweeper, since its
        /// already homogenized.
    }

    void PlaneSweeper_2D3D::calc_fission_source( float_t k, ArrayX
            &fission_source ) const {
        /// \todo implement. not sure what to do here (how should the
        /// EigenSolver see this sweeper?) Nice to have a choice for once, at
        /// least...
    }

    void PlaneSweeper_2D3D::output( H5File& file ) const {
        sn_sweeper_.output( file );
        moc_sweeper_.output( file );
    }
}
