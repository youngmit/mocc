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

#include "plane_sweeper_2d3d.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "util/range.hpp"

#include "core/error.hpp"

using mocc::sn::SnSweeper;

using std::cout;
using std::endl;
using std::cin;
using std::setfill;
using std::setw;

namespace mocc { namespace cmdo {
////////////////////////////////////////////////////////////////////////////////
    /// \todo make sure to check the angular quadratures for conformance
    PlaneSweeper_2D3D::PlaneSweeper_2D3D( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        TransportSweeper(input),
        mesh_(mesh),
        pair_(SnSweeperFactory_CDD(input.child("sn_sweeper"), mesh)),
        sn_sweeper_(std::move(pair_.first)),
        corrections_( pair_.second ),
        moc_sweeper_(input.child("moc_sweeper"), mesh),
        ang_quad_(moc_sweeper_.get_ang_quad()),
        tl_( sn_sweeper_->n_group(), mesh_.n_pin() ),
        sn_resid_( sn_sweeper_->n_group() ),
        prev_moc_flux_( sn_sweeper_->n_group(), mesh_.n_pin()),
        i_outer_( -1 )
    {
        this->parse_options( input );
        core_mesh_ = &mesh;

        xs_mesh_ = moc_sweeper_.get_xs_mesh();
        flux_.reference(moc_sweeper_.flux());
        vol_ = moc_sweeper_.volumes();

        n_reg_ = moc_sweeper_.n_reg();
        n_group_ = xs_mesh_->n_group();
        groups_ = Range(n_group_);

        auto sn_xs_mesh =
            sn_sweeper_->get_homogenized_xsmesh();
        assert( corrections_ );
        moc_sweeper_.set_coupling( corrections_, sn_xs_mesh );

        if( !keep_sn_quad_ ) {
            sn_sweeper_->set_ang_quad(ang_quad_);
        }

        sn_sweeper_->get_homogenized_xsmesh()->set_flux( moc_sweeper_.flux() );

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
        bool do_moc = ((i_outer_+1) > n_inactive_moc_) &&
            ((i_outer_ % moc_modulo_) == 0);
        if( do_moc ) {
            moc_sweeper_.sweep( group );

            int n_negative = 0;
            const auto flux = moc_sweeper_.flux()( blitz::Range::all(), group );
            for( const auto &v: flux ) {
                if( v < 0.0 ) {
                    n_negative++;
                }
            }
            if( n_negative > 0 ) {
                cout << n_negative << " negative fluxes in group " << group <<
                    endl;
            }
        }

        ArrayB1 prev_moc_flux = prev_moc_flux_( group,
                blitz::Range::all() );
        prev_moc_flux(0) = 0.0;
        moc_sweeper_.get_pin_flux_1g( group, prev_moc_flux );

        if( do_mocproject_ ) {
            sn_sweeper_->set_pin_flux_1g( group, prev_moc_flux );
        }

        // Sn sweeper
        // For now, we are assuming that the Sn cross sections are being updated
        // in the last inner iteration of the MoC. We should come up with
        // something better, since that actually ignores the change to the fine
        // mesh flux in the last MoC inner. For low numbers of MoC inners, this
        // essentially constitutes an outer-iteration lag of the updated cross
        // sections :-/
        sn_sweeper_->sweep( group );

        if( do_snproject_ ) {
            ArrayB1 sn_flux(mesh_.n_pin());
            sn_sweeper_->get_pin_flux_1g( group, sn_flux );
            moc_sweeper_.set_pin_flux_1g( group, sn_flux );
        }

        // Compute Sn-MoC residual
        real_t residual = 0.0;
        for( size_t i=0; i<prev_moc_flux.size(); i++ ) {
            residual += (prev_moc_flux(i) - sn_sweeper_->flux( group, i )) *
                        (prev_moc_flux(i) - sn_sweeper_->flux( group, i ));
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
        sn_sweeper_->initialize();
        moc_sweeper_.initialize();
    }

////////////////////////////////////////////////////////////////////////////////
    void PlaneSweeper_2D3D::get_pin_flux_1g( int ig, ArrayB1 &flux ) const {
        if( expose_sn_ ) {
            sn_sweeper_->get_pin_flux_1g( ig, flux );
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


        // Hand the transverse leakage to the MoC sweeper.
        moc_sweeper_.apply_transverse_leakage( group, tl_fsr );

    }

////////////////////////////////////////////////////////////////////////////////
    void PlaneSweeper_2D3D::output( H5Node &file ) const {
        // Put the Sn data in its own location
        {
            auto g = file.create_group( "/Sn" );
            sn_sweeper_->output( g );
        }

        // Put the MoC data in its own location
        {
            auto g = file.create_group( "/MoC" );
            moc_sweeper_.output( g );
        }

        VecI dims;
        dims.push_back(mesh_.nz());
        dims.push_back(mesh_.ny());
        dims.push_back(mesh_.nx());

        // Write out the Sn-MoC residual convergence
        file.create_group("/SnResid");
        for( int g=0; g<n_group_; g++ ) {
            std::stringstream setname;
            setname << "/SnResid/" << setfill('0') << setw(3) << g;
            VecI niter(1, sn_resid_[g].size());
            file.write( setname.str(), sn_resid_[g], niter );
        }

        {
            auto flux = prev_moc_flux_.copy();
            Normalize( flux.begin(), flux.end());
            auto h5g = file.create_group("moc_flux");
            for( const auto &ig: groups_ ) {
                std::stringstream setname;
                setname << setfill('0') << setw(3) << ig+1;
                h5g.write( setname.str(), flux(ig, blitz::Range::all()), dims );
            }
        }

        // Write out the transverse leakages
        {
            auto group = file.create_group("/transverse_leakage");
            for( int g=0; g<n_group_; g++ ) {
                std::stringstream setname;
                setname << setfill('0') << setw(3) << g;

                const auto tl_slice = tl_((int)g, blitz::Range::all());

                group.write( setname.str(), tl_slice.begin(), tl_slice.end(),
                        dims );
            }
        }

        // Write out the correction factors
        corrections_->output(file);
    }

////////////////////////////////////////////////////////////////////////////////
    // At some point it might be nice to make the options const and initialized
    // them in the initializer list, then just check for validity later. This if
    // fine for now.
    void PlaneSweeper_2D3D::parse_options( const pugi::xml_node &input ) {
        // Set defaults for everything
        expose_sn_ = false;
        do_snproject_ = false;
        do_mocproject_ = false;
        keep_sn_quad_ = false;
        do_tl_ = true;
        n_inactive_moc_ = 0;
        moc_modulo_ = 1;
        relax_ = 1.0;

        // Override with entries in the input node
        if( !input.attribute("expose_sn").empty() ) {
            expose_sn_ = input.attribute("expose_sn").as_bool();
        }
        if( !input.attribute("sn_project").empty() ) {
            do_snproject_ = input.attribute("sn_project").as_bool();
        }
        if( !input.attribute("moc_project").empty() ) {
            do_mocproject_ = input.attribute("moc_project").as_bool();
        }
        if( !input.attribute("tl").empty() ) {
            do_tl_ = input.attribute("tl").as_bool();
        }
        if( !input.attribute("inactive_moc").empty() ) {
            n_inactive_moc_ = input.attribute("inactive_moc").as_int();
        }
        if( !input.attribute("moc_modulo").empty() ) {
            moc_modulo_ = input.attribute("moc_modulo").as_int();
        }
        if( !input.attribute("preserve_sn_quadrature").empty() ) {
            keep_sn_quad_ = input.attribute("preserve_sn_quadrature").as_bool();
        }
        if( !input.attribute("relax").empty() ) {
            relax_ = input.attribute("relax").as_bool();
        }

        LogFile << "2D3D Sweeper options:" << std::endl;
        LogFile << "    Sn Projection: " << do_snproject_ << std::endl;
        LogFile << "    MoC Projection: " << do_mocproject_ << std::endl;
        LogFile << "    Expose Sn pin flux: " << expose_sn_ << std::endl;
        LogFile << "    Keep original Sn quadrature: " << keep_sn_quad_ << std::endl;
        LogFile << "    Transverse Leakage: " << do_tl_ << std::endl;
        LogFile << "    Inactive MoC Outer Iterations: "
            << n_inactive_moc_ << std::endl;
        LogFile << "    MoC sweep modulo: " << moc_modulo_ << std::endl;


    }
} } // Namespace mocc::cmdo
