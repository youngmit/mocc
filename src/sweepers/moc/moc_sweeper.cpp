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

#include "moc_sweeper.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>

#include "error.hpp"
#include "files.hpp"
#include "moc_current_worker.hpp"
#include "string_utils.hpp"
#include "utils.hpp"

using std::endl;
using std::cout;
using std::cin;

mocc::VecF temp;

/**
 * \brief Return the appropriate sizing values for construcing a \ref
 * mocc::BoundaryCondition.
 *
 * This exists so that it may be used to construct the boundary condition
 * members of the \ref mocc::moc::MoCSweeper from the initializer list.
 */
std::vector<mocc::BC_Size_t> bc_size_helper( const mocc::moc::RayData &rays )
{
    std::vector<mocc::BC_Size_t> bc_dims( rays.ang_quad().ndir_oct()*4 );
    for( int iang=0; iang<(int)rays.begin()->size(); iang++ ) {
        int iang1 = iang;
        int iang2 = rays.ang_quad().reverse(iang);
        int nx = rays.ny(iang);
        int ny = rays.nx(iang);
        bc_dims[iang1] = {nx, ny, 0};
        bc_dims[iang2] = {nx, ny, 0};
    }

    assert( bc_dims.size() == bc_dims.capacity() );
    return bc_dims;
}

namespace mocc { namespace moc {
    MoCSweeper::MoCSweeper( const pugi::xml_node& input,
            const CoreMesh& mesh ):
        TransportSweeper( input, mesh ),
        timer_(RootTimer.new_timer( "MoC Sweeper", true) ),
        timer_init_(timer_.new_timer( "Initialization", true) ),
        timer_sweep_(timer_.new_timer( "Sweep" )),
        mesh_( mesh ),
        rays_( input.child("rays"),
               ang_quad_,
               mesh
             ),
        boundary_( mesh.nz(),
            BoundaryCondition( n_group_,
                ang_quad_,
                mesh_.boundary(),
                bc_size_helper(rays_) ) ),
        boundary_out_( mesh.nz(),
            BoundaryCondition( 1,
                ang_quad_,
                mesh_.boundary(),
                bc_size_helper(rays_) ) ),
        xstr_( n_reg_ ),
        flux_1g_( ),
        bc_type_( mesh_.boundary() ),
        dump_rays_(false),
        gauss_seidel_boundary_(true),
        allow_splitting_(false)
    {

        LogFile << "Constructing a base MoC sweeper" << std::endl;

        // Make sure we have input from the XML
        if( input.empty() ) {
            throw EXCEPT("No input specified to initialize MoC sweeper.");
        }

        // Parse the number of inner iterations
        int int_in = input.attribute("n_inner").as_int(-1);
        if(int_in < 0) {
            throw EXCEPT("Invalid number of inner iterations specified "
                "(n_inner).");
        }
        n_inner_ = int_in;

        // Parse the output options
        dump_rays_ = input.attribute("dump_rays").as_bool(false);

        // Determine boundary update technique
        gauss_seidel_boundary_ = true;
        if( !input.attribute("boundary").empty() ) {
            std::string in_string = input.attribute("boundary").value();
            sanitize(in_string);
            if( in_string == "jacobi" || in_string == "j" ) {
                gauss_seidel_boundary_ = false;
            } else if ( in_string == "gs" ) {

            } else {
                throw EXCEPT("Unrecognized boundary update option.");
            }
        }

        // Parse TL source splitting setting
        allow_splitting_ = input.attribute("tl_splitting").as_bool(false);
        if( allow_splitting_ ) {
            split_.resize(n_reg_);
        }

        // Set up the array of volumes (surface area)
        int ireg = 0;
        for( auto &pin: mesh_ ) {
            const VecF& pin_vols = pin->vols();
            for( auto &v: pin_vols ) {
                vol_[ireg] = v;
                ireg++;
            }
        }

        if( dump_rays_ ) {
            std::ofstream rayfile("rays.py");
            rayfile << rays_ << endl;
        }

        // Replace the angular quadrature with the modularized version
        ang_quad_ = rays_.ang_quad();

        timer_init_.toc();
        timer_.toc();

        return;
    } // MoCSweeper( input, mesh )

    void MoCSweeper::sweep( int group ) {
        assert(source_);

        timer_.tic();
        timer_sweep_.tic();

        // set up the xstr_ array
        this->expand_xstr( group );

        flux_1g_.reference( flux_( blitz::Range::all(), group ) );

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            // update the self-scattering source
            source_->self_scatter( group, xstr_ );

            // Perform the stock sweep unless we are on the last outer and have
            // a CoarseData object.
            if( inner == n_inner_-1 && coarse_data_ ) {
                // Wipe out the existing currents (only on X- and Y-normal
                // faces)
                coarse_data_->zero_data_radial( group );

                moc::Current cw( coarse_data_, &mesh_ );
                this->sweep1g( group, cw );
                coarse_data_->set_has_radial_data(true);
            } else {
                moc::NoCurrent cw( coarse_data_, &mesh_ );
                this->sweep1g( group, cw );
            }
        }

        timer_.toc();
        timer_sweep_.toc();
        return;
    } // sweep( group )

    void MoCSweeper::initialize() {
        // Set the flux on the coarse mesh
        if( coarse_data_ ) {
            coarse_data_->flux = 1.0;
        }
        // There are better ways to do this, but for now, just start with 1.0
        flux_ = 1.0;

        // Walk through the boundary conditions and initialize them the 1/4pi
        real_t val = 1.0/FPI;
        for( auto &boundary: boundary_ ) {
            boundary.initialize_scalar(val);
        }

        return;
    } // initialize()

    void MoCSweeper::expand_xstr( int group ) {
        if( allow_splitting_ ) {
            for( auto &xsr: *xs_mesh_ ) {
                real_t xstr = xsr.xsmactr(group);
                for( auto &ireg: xsr.reg() ) {
                    xstr_(ireg) = xstr + split_(ireg);
                }
            }
        } else {
            for( auto &xsr: *xs_mesh_ ) {
                real_t xstr = xsr.xsmactr(group);
                for( auto &ireg: xsr.reg() ) {
                    xstr_(ireg) = xstr;
                }
            }
        }
    } // expand_xstr( group )

    void MoCSweeper::get_pin_flux_1g( int group, ArrayB1 &flux ) const {
        assert(flux.size() == mesh_.n_pin() );
        for( auto &f: flux ) {
            f = 0.0;
        }

        int ireg = 0;
        int ipin = 0;
        for( auto &pin: mesh_ ) {
            Position pos = mesh_.pin_position(ipin);
            int i = mesh_.coarse_cell(pos);
            real_t v = 0.0;
            for( int ir=0; ir<pin->n_reg(); ir++ ) {
                v += vol_[ireg];
                flux(i) += flux_(ireg, group)*vol_[ireg];
                ireg++;
            }
            flux(i) /= v;
            ipin++;
        }

        return;
    }

    real_t MoCSweeper::set_pin_flux_1g( int group,
            const ArrayB1 &pin_flux)
    {
        real_t resid = 0.0;
        size_t ipin = 0;
        int ireg = 0;
        for( const auto &pin: mesh_ ) {
            size_t i_coarse = mesh_.coarse_cell( mesh_.pin_position(ipin) );
            real_t fm_flux = 0.0;

            int stop = ireg+pin->n_reg();
            for( ; ireg<stop; ireg++ ) {
                fm_flux += vol_[ireg]*flux_(ireg, group);
            }

            fm_flux /= pin->vol();
            real_t f = pin_flux(i_coarse)/fm_flux;

            ireg -= pin->n_reg();
            for( ; ireg<stop; ireg++ ) {
                flux_(ireg, group) = flux_(ireg, group) * f;
            }

            real_t e = fm_flux - pin_flux(i_coarse);
            resid += e*e;
            ipin++;
        }
        return std::sqrt(resid);
    } // set_pin_flux_1f( group, pin_flux )

    void MoCSweeper::apply_transverse_leakage( int group, const ArrayB1 &tl ) {
        assert((int)tl.size() == n_reg_);

        flux_1g_.reference( flux_( blitz::Range::all(), group ) );
        
        /// \todo for now, this is using a pretty invasive direct access the the
        /// source. Might be good to do as a call to auxiliary() instead
        if( allow_splitting_ ) {
            int n_split = 0;
            split_ = 0.0;
            for( int ireg = 0; ireg < n_reg_; ireg++ ) {
                real_t s = (*source_)[ireg] + tl(ireg);
                /// \todo once we get this working, clean this branch up
                if( s < 0.0 ) {
                    n_split++;
                    split_(ireg) = -s/flux_1g_(ireg); 
                    (*source_)[ireg] = 0.0;
                } else {
                    (*source_)[ireg] = s;
                }
            }

            if( n_split > 0 ) {
                LogScreen << "Split " << n_split << " region sources"
                          << std::endl;
            }
        } else {
            for( int ireg = 0; ireg < n_reg_; ireg++ ) {
                real_t s = (*source_)[ireg] + tl(ireg);
                (*source_)[ireg] = s;
            }
        }

        return;
    } // apply_transverse_leakage( tl )

    /**
     * \brief Check for the balance of neutrons within each pin cell.
     *
     * \todo Make sure this is valid in the presence of source splitting
     */
    void MoCSweeper::check_balance( int group ) const {
        ArrayB1 b(mesh_.n_pin());
        b = 0.0;

        // Get the removal cross section in a nice format
        ArrayB1 xsrm(n_reg_);
        xsrm = 0.0;
        for( auto &xsr: *xs_mesh_ ) {
            real_t rm = xsr.xsmacrm(group);
            for( auto &ireg: xsr.reg() ) {
                xsrm(ireg) = rm;
            }
        }

        ArrayB1 current_1g = coarse_data_->current(blitz::Range::all(), group);

        int ipin = 0;
        int ireg = 0;
        for( const auto pin: mesh_ ) {
            int icell = mesh_.coarse_cell(mesh_.pin_position(ipin));
            real_t bi = 0.0;

            for( int ireg_pin=0; ireg_pin<pin->n_reg(); ireg_pin++ ) {
                bi -= flux_(ireg, group)*vol_[ireg]*xsrm(ireg);
                bi += (*source_)[ireg]*vol_[ireg];
                ireg++;
            }

            // Current
            bi -= current_1g( mesh_.coarse_surf(icell, Surface::EAST) ) *
                mesh_.coarse_area(icell, Surface::EAST);
            bi -= current_1g( mesh_.coarse_surf(icell, Surface::NORTH) ) *
                mesh_.coarse_area(icell, Surface::NORTH);
            bi -= current_1g( mesh_.coarse_surf(icell, Surface::TOP) ) *
                mesh_.coarse_area(icell, Surface::TOP);
            bi += current_1g( mesh_.coarse_surf(icell, Surface::WEST) ) *
                mesh_.coarse_area(icell, Surface::WEST);
            bi += current_1g( mesh_.coarse_surf(icell, Surface::SOUTH) ) *
                mesh_.coarse_area(icell, Surface::SOUTH);
            bi += current_1g( mesh_.coarse_surf(icell, Surface::BOTTOM) ) *
                mesh_.coarse_area(icell, Surface::BOTTOM);

            b(icell) = bi;
            ipin++;
        }

        cout << "MoC cell balance:" << endl;
        for( auto v: b ) {
            cout << v << endl;
        }

        return;
    }

    void MoCSweeper::output( H5Node &node ) const {
        // Get core dimensions from the mesh
        VecI dims = mesh_.dimensions();
        std::reverse( dims.begin(), dims.end() );

        // Make a group in the file to store the flux
        node.create_group("flux");

        ArrayB2 flux = this->get_pin_flux();
        Normalize( flux.begin(), flux.end() );

        LogFile << "Boundary update: ";
        if( gauss_seidel_boundary_ ) {
            LogFile << "Gauss-Seidel" << std::endl;
        } else {
            LogFile << "Jacobi" << std::endl;
        }

        for( int ig=0; ig<n_group_; ig++ ){
            std::stringstream setname;
            setname << "flux/" << std::setfill('0') << std::setw(3) << ig+1;

            ArrayB1 flux_1g = flux(blitz::Range::all(), ig);
            node.write( setname.str(), flux_1g.begin(),
                    flux_1g.end(), dims );
        }

        // Pin powers
        auto pow = this->pin_powers();
        node.write("pin_powers", pow);

        ang_quad_.output( node );

        return;
    }
} }
