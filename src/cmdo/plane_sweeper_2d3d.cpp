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
        mesh_(mesh),
        sn_sweeper_(input.child("sn_sweeper"), mesh),
        moc_sweeper_(input.child("moc_sweeper"), mesh),
        ang_quad_(moc_sweeper_.get_ang_quad())
    {
        /// \todo initialize the rest of the data members on the TS base type
        xs_mesh_ = moc_sweeper_.get_xs_mesh();
        n_reg_ = moc_sweeper_.n_reg();
        n_group_ = xs_mesh_->n_group();

        corrections_ = CorrectionData( sn_sweeper_.n_reg(), ang_quad_.ndir()/2,
                n_group_ );

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

        sn_sweeper_.get_homogenized_xsmesh()->update( moc_sweeper_.flux() );
        sn_sweeper_.sweep( group );

        // Compute Sn-MoC residual
        VecF moc_flux;
        moc_sweeper_.get_pin_flux_1g( group, moc_flux );

        real_t residual = 0.0;
        for( size_t i=0; i<moc_flux.size(); i++ ) {
            residual += (moc_flux[i] - sn_sweeper_.flux( group, i )) *
                        (moc_flux[i] - sn_sweeper_.flux( group, i ));
        }
        residual = sqrt(residual);
cout << "MoC/Sn residual: " << residual << endl;
    }

    void PlaneSweeper_2D3D::initialize() {
        sn_sweeper_.initialize();
        moc_sweeper_.initialize();
    }

    void PlaneSweeper_2D3D::get_pin_flux_1g( int ig, VecF &flux ) const {
        sn_sweeper_.get_pin_flux_1g( ig, flux );
    }

    VecF PlaneSweeper_2D3D::get_pin_flux() const {
        VecF flux( mesh_.n_pin(), 0.0 );

        assert(false);

        return flux;
    }

    real_t PlaneSweeper_2D3D::total_fission( bool old ) const {
        // using the MoC for now, but once things start to converge, might want
        // to use Sn, since it should yield the same result at convergence and
        // is cheaper to evaluate.
        real_t tfis = moc_sweeper_.total_fission( old );
        return tfis;
    }

    void PlaneSweeper_2D3D::calc_fission_source( real_t k, 
            ArrayF &fission_source ) const {
        moc_sweeper_.calc_fission_source( k, fission_source );
        return;
    }

    void PlaneSweeper_2D3D::output( H5::CommonFG *file ) const {
        // Put the Sn data in its own location
        {
            H5::Group g = file->createGroup( "/Sn" );
            sn_sweeper_.output( &g );
        }

        // Put the MoC data in its own location
        {
            H5::Group g = file->createGroup( "/MoC" );
            moc_sweeper_.output( &g );
        }

        // Write out the correction factors
        file->createGroup("/alpha_x");
        file->createGroup("/alpha_y");
        file->createGroup("/beta");

        VecI dims;
        dims.push_back(mesh_.nz());
        dims.push_back(mesh_.ny());
        dims.push_back(mesh_.nx());
        for( size_t g=0; g<n_group_; g++ ) {
            for( size_t a=0; a<ang_quad_.ndir_oct()*4; a++ ) {
                VecF alpha_x( mesh_.n_pin(), 0.0 );
                VecF alpha_y( mesh_.n_pin(), 0.0 );
                VecF beta( mesh_.n_pin(), 0.0 );
                for( size_t i=0; i<mesh_.n_pin(); i++ ) {
                    alpha_x[i] = corrections_.alpha( i, a, g, Normal::X_NORM );
                    alpha_y[i] = corrections_.alpha( i, a, g, Normal::Y_NORM );
                    beta[i] = corrections_.beta( i, a, g );
                }
                
                {
                    std::stringstream setname;
                    setname << "/beta/" << g << "_" << a;
                    HDF::Write( file, setname.str(), beta, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_x/" << g << "_" << a;
                    HDF::Write( file, setname.str(), alpha_x, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_y/" << g << "_" << a;
                    HDF::Write( file, setname.str(), alpha_y, dims );
                }
            }
        }
    }
}
