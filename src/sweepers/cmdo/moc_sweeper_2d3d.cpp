#include "moc_sweeper_2d3d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <valarray>

#include "correction_worker.hpp"
#include "files.hpp"

using std::cout;
using std::cin;
using std::endl;

namespace mocc { namespace cmdo {
    MoCSweeper_2D3D::MoCSweeper_2D3D( const pugi::xml_node &input,
            const CoreMesh &mesh ):
        MoCSweeper( input, mesh ),
        corrections_( nullptr ),
        sn_xs_mesh_( nullptr ),
        internal_coupling_( false ),
        correction_residuals_( n_group_ )
    {
        LogFile << "Constructing a 2D3D MoC sweeper" << std::endl;
        return;
    };

    void MoCSweeper_2D3D::sweep( int group ) {
        timer_.tic();
        timer_sweep_.tic();
        assert(source_);
        assert(sn_xs_mesh_);

        if( !coarse_data_ ) {
            throw EXCEPT("2D3D MoC sweeper needs coarse data to collect "
                    "calculate correction factors. Try enabling CMFD.");
        }

        n_sweep_++;

        // set up the xstr_ array
        for( auto &xsr: *xs_mesh_ ) {
            real_t xstr = xsr.xsmactr()[group];
            for( auto &ireg: xsr.reg() ) {
                xstr_(ireg) = xstr;
            }
        }

        // Instantiate the workers for current/no current
        CurrentCorrections ccw( coarse_data_, &mesh_, corrections_.get(),
                source_->get_transport( 0 ), xstr_, ang_quad_, *sn_xs_mesh_,
                rays_ );
        moc::NoCurrent ncw( coarse_data_, &mesh_ );


        flux_1g_.reference( flux_( blitz::Range::all(), group ) );

        // Perform inner iterations
        for( unsigned int inner=0; inner<n_inner_; inner++ ) {
            n_sweep_inner_++;
            // update the self-scattering source
            source_->self_scatter( group );

            // Perform the stock sweep unless we are on the last outer and have
            // a CoarseData object.
            if( inner == n_inner_-1 && coarse_data_ ) {
                coarse_data_->zero_data_radial( group );
                sn_xs_mesh_->update();
                this->sweep1g( group, ccw );
                coarse_data_->set_has_radial_data(true);
                correction_residuals_[group].push_back(ccw.residual());
            } else {
                this->sweep1g( group, ncw );
            }
        }

        timer_.toc();
        timer_sweep_.toc();
        return;
    }

    void MoCSweeper_2D3D::output( H5Node &node ) const {
        LogFile << "MoC Sweeper 2D3D output:" << std::endl;
        LogFile << "    Number of sweeps, outer: " << n_sweep_ << std::endl;
        LogFile << "    Number of sweeps, inner: " << n_sweep_inner_
                << std::endl;;

        MoCSweeper::output( node );
        if( internal_coupling_ ) {
            corrections_->output( node );
            sn_xs_mesh_->update();
            sn_xs_mesh_->output( node );
        }
        
        auto residual_group = node.create_group("correction_residual");

        for( int ig=0; ig<n_group_; ig++ ) {
            std::stringstream g_str;
            g_str << ig+1;
            auto g = residual_group.create_group(g_str.str());
            VecF ax;
            VecF ay;
            VecF b;
            for( auto data: correction_residuals_[ig] ) {
                ax.push_back(data[0]);
                ay.push_back(data[1]);
                b.push_back(data[2]);
            }

            g.write("alpha_x", ax);
            g.write("alpha_y", ay);
            g.write("beta", b);
        }
    }
} } // Namespace mocc::cmdo
