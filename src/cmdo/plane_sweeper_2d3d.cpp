#include "plane_sweeper_2d3d.hpp"

#include <cmath>

#include "error.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
    /// \todo make sure to check the angular quadratures for conformance
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
        const XSMeshHomogenized* sn_xs_mesh = 
            sn_sweeper_.get_homogenized_xsmesh().get();
        moc_sweeper_.set_coupling( &corrections_, sn_xs_mesh );

        coarse_data_ = nullptr;
        
        return;
    }

    void PlaneSweeper_2D3D::sweep( int group ) {
        if (!coarse_data_) {
            throw EXCEPT("CMFD must be enabled to do 2D3D.");
        }
        moc_sweeper_.sweep( group );
        sn_sweeper_.sweep( group );

        // Compute Sn-MoC residual
        VecF moc_flux;
        moc_sweeper_.get_pin_flux( group, moc_flux );
cout << "representative moc flux: " << moc_flux[0] << endl;
        float_t residual = 0.0;
        for( size_t i=0; i<moc_flux.size(); i++ ) {
            cout << moc_flux[i] << " " << sn_sweeper_.flux( group, i ) << endl;
            residual += (moc_flux[i] - sn_sweeper_.flux( group, i )) *
                        (moc_flux[i] - sn_sweeper_.flux( group, i ));
        }
        residual = sqrt(residual);
        cout << "MoC/Sn residual: " << residual << endl;
        cin.ignore();
        
    }

    void PlaneSweeper_2D3D::initialize() {
        sn_sweeper_.initialize();
        moc_sweeper_.initialize();
    }

    void PlaneSweeper_2D3D::get_pin_flux( int ig, VecF &flux ) const {
        sn_sweeper_.get_pin_flux( ig, flux );
    }

    float_t PlaneSweeper_2D3D::total_fission( bool old ) const {
        // using the MoC for now, but once things start to converge, might want
        // to use Sn, since it should yield the same result at convergence and
        // is cheaper to evaluate.
        float_t tfis = moc_sweeper_.total_fission( old );
        return tfis;
    }

    void PlaneSweeper_2D3D::calc_fission_source( float_t k, 
            ArrayX &fission_source ) const {
        moc_sweeper_.calc_fission_source( k, fission_source );
        return;
    }

    void PlaneSweeper_2D3D::output( H5File& file ) const {
        // We need to be a little careful about how we do output to avoid
        // collisions in the HDF5 tree. For now, lets just output Sn data.
        sn_sweeper_.output( file );
        //moc_sweeper_.output( file );
    }
}