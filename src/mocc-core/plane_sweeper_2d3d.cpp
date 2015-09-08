#include "plane_sweeper_2d3d.hpp"

namespace mocc {
    PlaneSweeper_2D3D::PlaneSweeper_2D3D( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        sn_sweeper_(input.child("sn_sweeper"), mesh),
        moc_sweeper_(input.child("moc_sweeper"), mesh)
    {
        /// \todo initialize the rest of the data members on the TS base type
        return;
    }

    void PlaneSweeper_2D3D::sweep( int group ) {
        /// \todo implement
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
        /// \todo implement.
    }
}
