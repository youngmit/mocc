#include "plane_sweeper_2d3d.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "error.hpp"

using mocc::sn::SnSweeper;

using std::cout;
using std::endl;
using std::cin;
using std::setfill;
using std::setw;

namespace mocc {
////////////////////////////////////////////////////////////////////////////////
    /// \todo make sure to check the angular quadratures for conformance
    PlaneSweeper_2D3D::PlaneSweeper_2D3D( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        TransportSweeper(input),
        mesh_(mesh),
        sn_sweeper_(input.child("sn_sweeper"), mesh),
        moc_sweeper_(input.child("moc_sweeper"), mesh),
        ang_quad_(moc_sweeper_.get_ang_quad()),
        corrections_( new CorrectionData( sn_sweeper_.n_reg(),
            ang_quad_.ndir()/2,
            sn_sweeper_.n_group() ) ),
        tl_( sn_sweeper_.n_group(), mesh_.n_pin() ),
        sn_resid_( sn_sweeper_.n_group() ),
        i_outer_( 0 )
    {
        this->parse_options( input );
        /// \todo initialize the rest of the data members on the TS base type
        core_mesh_ = &mesh;

        if( expose_sn_ ) {
            xs_mesh_ = sn_sweeper_.get_xs_mesh();
            n_reg_ = sn_sweeper_.n_reg();
            flux_.reference(sn_sweeper_.flux());
            vol_ = sn_sweeper_.volumes();
        } else {
            xs_mesh_ = moc_sweeper_.get_xs_mesh();
            n_reg_ = moc_sweeper_.n_reg();
            flux_.reference(moc_sweeper_.flux());
            vol_ = moc_sweeper_.volumes();
        }

        n_group_ = xs_mesh_->n_group();

        sn_sweeper_.worker()->set_corrections( corrections_ );
        const XSMeshHomogenized* sn_xs_mesh =
            sn_sweeper_.get_homogenized_xsmesh().get();
        moc_sweeper_.set_coupling( corrections_, sn_xs_mesh );

        sn_sweeper_.set_ang_quad(ang_quad_);

        sn_sweeper_.get_homogenized_xsmesh()->set_flux( flux_ );

        coarse_data_ = nullptr;


        return;
    }

////////////////////////////////////////////////////////////////////////////////
    void PlaneSweeper_2D3D::sweep( int group ) {
        if( !coarse_data_ ) {
            throw EXCEPT("CMFD must be enabled to do 2D3D.");
        }

        /// \todo do something less brittle
        if( group == 0 ) {
            i_outer_++;
        }

        // Calculate transverse leakage source
        if( do_tl_ ) {
            this->add_tl( group );
        }

        // MoC Sweeper
        if( i_outer_ > n_inactive_moc_ ) {
            moc_sweeper_.sweep( group );
        }
        ArrayB1 moc_flux( mesh_.n_pin() );
        moc_sweeper_.get_pin_flux_1g( group, moc_flux );

        // Sn sweeper
        sn_sweeper_.get_homogenized_xsmesh()->update( );
        sn_sweeper_.sweep( group );

        if( do_snproject_ ) {
            ArrayB1 sn_flux(mesh_.n_pin());
            sn_sweeper_.get_pin_flux_1g( group, sn_flux );
            moc_sweeper_.set_pin_flux_1g( group, sn_flux );
        }

        // Compute Sn-MoC residual

        real_t residual = 0.0;
        for( size_t i=0; i<moc_flux.size(); i++ ) {
            residual += (moc_flux(i) - sn_sweeper_.flux( group, i )) *
                        (moc_flux(i) - sn_sweeper_.flux( group, i ));
        }
        residual = sqrt(residual)/mesh_.n_pin();

        cout << "MoC/Sn residual: " << residual;
        if( sn_resid_[group].size() > 0 ) {
             cout << "   \t" << residual-sn_resid_[group].back();
        }
        cout << endl;
        sn_resid_[group].push_back(residual);
    }

////////////////////////////////////////////////////////////////////////////////
    void PlaneSweeper_2D3D::initialize() {
        sn_sweeper_.initialize();
        moc_sweeper_.initialize();
    }

////////////////////////////////////////////////////////////////////////////////
    void PlaneSweeper_2D3D::get_pin_flux_1g( int ig, ArrayB1 &flux ) const {
        if( expose_sn_ ) {
            sn_sweeper_.get_pin_flux_1g( ig, flux );
        } else {
            moc_sweeper_.get_pin_flux_1g( ig, flux );
        }
    }

////////////////////////////////////////////////////////////////////////////////
    void PlaneSweeper_2D3D::add_tl( int group ) {
        assert( coarse_data_ );
        ArrayB1 tl_fsr( n_reg_ );

        int ireg_pin = 0;
        int ipin = 0;

        blitz::Array<real_t, 1> tl_g = tl_(group, blitz::Range::all());

        for( const auto &pin: mesh_ ) {
            Position pos = mesh_.pin_position( ipin );
            size_t icell = mesh_.coarse_cell( pos );
            real_t dz = mesh_.dz( pos.z );
            int surf_up = mesh_.coarse_surf( icell, Surface::TOP );
            int surf_down = mesh_.coarse_surf( icell, Surface::BOTTOM );
            real_t j_up = coarse_data_->current(surf_up, group);
            real_t j_down = coarse_data_->current(surf_down, group);
            tl_g(ipin) = ( j_down - j_up ) / dz;

            for( int ir=0; ir<pin->n_reg(); ir++ ) {
                tl_fsr( ir+ireg_pin ) = tl_g(ipin);
            }
            ipin++;
            ireg_pin += pin->n_reg();
        }

        // Can add the TL as an auxiliary source directly to the Source_2D3D,
        // since it extends the MoC source in the first place
        source_->auxiliary( tl_fsr );
    }

////////////////////////////////////////////////////////////////////////////////
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

        VecI dims;
        dims.push_back(mesh_.nz());
        dims.push_back(mesh_.ny());
        dims.push_back(mesh_.nx());

        // Write out the Sn-MoC residual convergence
        file->createGroup("/SnResid");
        for( size_t g=0; g<n_group_; g++ ) {
            std::stringstream setname;
            setname << "/SnResid/" << setfill('0') << setw(3) << g;
            VecI niter(1, sn_resid_[g].size());
            HDF::Write( file, setname.str(), sn_resid_[g], niter );
        }

        // Write out the transverse leakages
        file->createGroup("/TL");
        for( size_t g=0; g<n_group_; g++ ) {
            std::stringstream setname;
            setname << "/TL/" << setfill('0') << setw(3) << g;

            const auto tl_slice = tl_((int)g, blitz::Range::all());

            HDF::Write( file, setname.str(), tl_slice.begin(), tl_slice.end(),
                    dims );
        }

        // Write out the correction factors
        file->createGroup("/alpha_x");
        file->createGroup("/alpha_y");
        file->createGroup("/beta");

        for( size_t g=0; g<n_group_; g++ ) {
            for( int a=0; a<ang_quad_.ndir_oct()*4; a++ ) {
                VecF alpha_x( mesh_.n_pin(), 0.0 );
                VecF alpha_y( mesh_.n_pin(), 0.0 );
                VecF beta( mesh_.n_pin(), 0.0 );
                for( size_t i=0; i<mesh_.n_pin(); i++ ) {
                    alpha_x[i] = corrections_->alpha( i, a, g, Normal::X_NORM );
                    alpha_y[i] = corrections_->alpha( i, a, g, Normal::Y_NORM );
                    beta[i] = corrections_->beta( i, a, g );
                }

                {
                    std::stringstream setname;
                    setname << "/beta/" << g << "_" << a;
                    HDF::Write( file, setname.str(), beta, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_x/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                    HDF::Write( file, setname.str(), alpha_x, dims );
                }

                {
                    std::stringstream setname;
                    setname << "/alpha_y/" << setfill('0') << setw(3) << g
                            << "_"         << setfill('0') << setw(3) << a;
                    HDF::Write( file, setname.str(), alpha_y, dims );
                }
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////
    // At some point it might be nice to make the options const and initialized
    // them in the initializer list, then just check for validity later. This if
    // fine for now.
    void PlaneSweeper_2D3D::parse_options( const pugi::xml_node &input ) {
        // Set defaults for everything
        expose_sn_ = false;
        do_snproject_ = true;
        do_tl_ = true;
        n_inactive_moc_ = 0;

        // Override with entries in the input node
        if( !input.attribute("expose_sn").empty() ) {
            expose_sn_ = input.attribute("expose_sn").as_bool();
        }
        if( !input.attribute("sn_project").empty() ) {
            do_snproject_ = input.attribute("sn_project").as_bool();
        }
        if( !input.attribute("tl").empty() ) {
            do_tl_ = input.attribute("tl").as_bool();
        }
        if( !input.attribute("inactive_moc").empty() ) {
            n_inactive_moc_ = input.attribute("inactive_moc").as_int();
        }

        LogFile << "2D3D Sweeper options:" << std::endl;
        LogFile << "    Sn Projection: " << do_snproject_ << std::endl;
        LogFile << "    Expose Sn variables: " << expose_sn_ << std::endl;
        LogFile << "    Transversd Leakage: " << do_snproject_ << std::endl;
        LogFile << "    Inactive MoC Outer Iterations: "
            << n_inactive_moc_ << std::endl;


    }
}
