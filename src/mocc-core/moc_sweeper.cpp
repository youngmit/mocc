#include "moc_sweeper.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>

#include "error.hpp"
#include "files.hpp"
#include "moc_current_worker.hpp"
#include "utils.hpp"

using std::endl;
using std::cout;
using std::cin;

mocc::VecF temp;

namespace mocc {
    MoCSweeper::MoCSweeper( const pugi::xml_node& input,
            const CoreMesh& mesh ):
        TransportSweeper( input, mesh ),
        mesh_( mesh ),
        rays_( input.child("rays"),
               ang_quad_,
               mesh
             ),
        xstr_( n_reg_ ),
        flux_1g_( n_reg_ ),
        qbar_( n_reg_ ),
        bc_type_( mesh_.boundary() )
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

        // Set up the array of volumes (surface area)
        int ireg = 0;
        for( auto &pin: mesh_ ) {
            const VecF& pin_vols = pin->vols();
            for( auto &v: pin_vols ) {
                vol_[ireg] = v;
                ireg++;
            }
        }

        // allocate space to store the boundary conditions
        boundary_.resize( n_group_ );
        for( auto &group_rays: boundary_ ) {
            group_rays.resize( mesh_.nz() );
            for( auto &angle_rays: group_rays ) {
                // We actually allocate BCs for all 4 octants to make things a
                // little simpler.
                angle_rays.resize( ang_quad_.ndir_oct()*4 );
                int iang = 0;
                for( auto ang_it=ang_quad_.octant(1);
                        ang_it!=ang_quad_.octant(5); ang_it++ ) {
                    angle_rays[iang].resize(rays_.n_rays(iang));
                    iang++;
                }
            }
        }
        boundary_out_.resize( mesh_.nz() );
        for( auto &angle_rays: boundary_out_ ) {
            // We actually allocate BCs for all 4 octants to make things a
            // little simpler.
            angle_rays.resize( ang_quad_.ndir_oct()*4 );
            int iang = 0;
            for( auto ang_it=ang_quad_.octant(1);
                    ang_it!=ang_quad_.octant(5); ang_it++ ) {
                angle_rays[iang].resize(rays_.n_rays(iang));
                iang++;
            }
        }

        if( dump_rays_ ) {
            std::ofstream rayfile("rays.py");
            rayfile << rays_ << endl;
        }

        return;
    }

    void MoCSweeper::sweep( int group ) {
        assert(source_);

        // set up the xstr_ array
        for( auto &xsr: *xs_mesh_ ) {
            real_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_[ireg] = xstr;
            }
        }

        flux_1g_ = flux_( blitz::Range::all(), group );

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            // update the self-scattering source
            source_->self_scatter( group, flux_1g_, qbar_ );

            // Perform the stock sweep unless we are on the last outer and have
            // a CoarseData object.
            if( inner == n_inner_-1 && coarse_data_ ) {
                // Wipe out the existing currents (only on X- and Y-normal
                // faces)
                this->zero_current( group ); 

                moc::Current cw( coarse_data_, &mesh_ );
                this->sweep1g( group, cw );
                coarse_data_->set_has_radial_data(true);
            } else {
                moc::NoCurrent cw( coarse_data_, &mesh_ );
                this->sweep1g( group, cw );
            }
        }

        flux_( blitz::Range::all(), group ) = flux_1g_;

        return;
    }


    void MoCSweeper::update_boundary( int group ) {
        int iplane=0;
        for( auto &plane_bcs: boundary_[group] ) {
            for( unsigned int iang=0; iang<plane_bcs.size(); iang++ ) {
                int nx = rays_.nx(iang);
                int ny = rays_.ny(iang);

                // Determine based on the angle quadrant which surfaces are
                // upwind and need to be updated for the given angle
                Surface upwind[2];
                if(ang_quad_[iang].ox > 0.0) {
                    upwind[0] = Surface::WEST;
                } else {
                    upwind[0] = Surface::EAST;
                }
                if(ang_quad_[iang].oy > 0.0) {
                    upwind[1] = Surface::SOUTH;
                } else {
                    upwind[1] = Surface::NORTH;
                }

                // Treat the two surfaces that were determined above
                if( bc_type_[(int)upwind[0]] == Boundary::REFLECT ) {
                    int ang_ref = ang_quad_.reflect(iang, upwind[0]);
                    for( int ibc=0; ibc<ny; ibc++ ) {
                        plane_bcs[iang][ibc] =
                            boundary_out_[iplane][ang_ref][ibc];
                    }
                } else {
                    for( int ibc=0; ibc<ny; ibc++ ) {
                        plane_bcs[iang][ibc] = 0.0;
                    }
                }

                if( bc_type_[(int)upwind[1]] == Boundary::REFLECT ) {
                    int ang_ref = ang_quad_.reflect(iang, upwind[1]);
                    for( int ibc=ny; ibc<(nx+ny); ibc++ ) {
                        plane_bcs[iang][ibc] =
                            boundary_out_[iplane][ang_ref][ibc];
                    }
                } else {
                    for( int ibc=ny; ibc<(nx+ny); ibc++ ) {
                        plane_bcs[iang][ibc] = 0.0;
                    }
                }
            } // angle loop
            iplane++;
        } // plane loop
    }

    void MoCSweeper::initialize() {
        // Set the flux on the coarse mesh
        if( coarse_data_ ) {
            coarse_data_->flux = 1.0;
        }
        // There are better ways to do this, but for now, just start with 1.0
        flux_ = 1.0;

        // Walk through the boundary conditions and initialize them the 1/4pi
        real_t val = 1.0/FPI;
        for( auto &group_rays: boundary_ ) {
            for( auto &plane_rays: group_rays ) {
                for( auto &angle_rays: plane_rays ) {
                    for( auto &ray: angle_rays ) {
                        ray = val;
                    }
                }
            }
        }
        return;
    }

    void MoCSweeper::get_pin_flux_1g( int group, ArrayB1 &flux ) const {
        assert(flux.size() == mesh_.n_pin() );
        for( auto &f: flux ) {
            f = 0.0;
        }

        int ireg = 0;
        int ipin = 0;
        for( auto &pin: mesh_ ) {
            Position pos = mesh_.pin_position(ipin);
            int i = mesh_.index_lex(pos);
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
            size_t i_coarse = mesh_.index_lex( mesh_.pin_position(ipin) );
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
    }

    void MoCSweeper::check_balance( int group ) const {
        ArrayF b(0.0, mesh_.n_pin());

        // Get the removal cross section in a nice format
        ArrayF xsrm(0.0, n_reg_);
        for( auto &xsr: *xs_mesh_ ) {
            real_t rm = xsr.xsmacrm()[group];
            for( auto &ireg: xsr.reg() ) {
                xsrm[ireg] = rm;
            }
        }

        ArrayB1 current_1g = coarse_data_->current(blitz::Range::all(), group);

        int ipin = 0;
        int ireg = 0;
        for( const auto pin: mesh_ ) {
            int icell = mesh_.coarse_cell(mesh_.pin_position(ipin));
            real_t bi = 0.0;

            for( int ireg_pin=0; ireg_pin<pin->n_reg(); ireg_pin++ ) {
                bi -= flux_(ireg, group)*vol_[ireg]*xsrm[ireg];
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

            b[icell] = bi;
            ipin++;
        }

        cout << "MoC cell balance:" << endl;
        for( auto v: b ) {
            cout << v << endl;
        }

        return;
    }

    void MoCSweeper::output( H5::CommonFG *node ) const {
        // Get core dimensions from the mesh
        VecI dims = mesh_.dimensions();
        std::reverse( dims.begin(), dims.end() );

        // Make a group in the file to store the flux
        node->createGroup("flux");

        ArrayB2 flux = this->get_pin_flux();
        Normalize( flux.begin(), flux.end() );

        for( size_t ig=0; ig<n_group_; ig++ ){
            std::stringstream setname;
            setname << "flux/" << std::setfill('0') << std::setw(3) << ig+1;

            ArrayB1 flux_1g = flux(blitz::Range::all(), ig);
            HDF::Write( node, setname.str(), flux_1g.begin(),
                    flux_1g.end(), dims );
        }

        return;
    }

    void MoCSweeper::zero_current( int group ) {
        assert( coarse_data_ );
        ArrayB1 current =
            coarse_data_->current(blitz::Range::all(), group);
        for( size_t plane=0; plane<mesh_.nz(); plane++ ) {
            for( auto surf=mesh_.plane_surf_xy_begin(plane);
                    surf!=mesh_.plane_surf_end(plane);
                    ++surf )
            {
                current(surf) = 0.0;
            }
        }
        return;
    }
}
