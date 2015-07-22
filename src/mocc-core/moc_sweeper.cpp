#include "moc_sweeper.hpp"

#include "error.hpp"

using std::endl;
using std::cout;

namespace mocc {
    
    MoCSweeper::MoCSweeper( const pugi::xml_node &input, 
            const CoreMesh &mesh ):
        TransportSweeper( mesh ),
        ang_quad_( input.child("ang_quad") ),
        rays_( input.child("rays"), ang_quad_, mesh )
    {   
        // Make sure we have input from the XML
        if( input.empty() ) {
            Error("No input specified to initialize MoC sweeper.");
        }

        // allocate space to store the boudary conditions
        boundary_.resize( ng_ );
        for( auto &group_rays: boundary_ ) {
            group_rays.resize( mesh_.nz() );
            for( auto &angle_rays: group_rays ) {
                angle_rays.resize( ang_quad_.ndir_oct()*2 );
                int iang = 0;
                for( auto ang_it=ang_quad_.octant(1); 
                        ang_it!=ang_quad_.octant(3); ang_it++ ) {
                    angle_rays[iang].resize(rays_.n_rays(iang));
                    iang++;
                }
            }
        }

        return;
    }

    void MoCSweeper::sweep( int group ) {
        int iplane = 0;
        for( auto &plane_rays: rays_ ) {
            int iang = 0;
            for( auto &ang_rays: plane_rays ) {
                int iray = 0;
                for( auto &ray: ang_rays ) {
                    // Initialize from bc
                    float_t psi = boundary_[group][iplane][iang][iray].fw;

                    // Propagate through core geometry
                    for( int iseg=0; iseg<ray.nseg(); iseg++ ) {
                        
                    }
                    
                    // Store boundary condition
                    
                    iray++;
                } // Rays
                iang++;
            } // angles
            iplane++;
        } // planes
        return;
    }

    void MoCSweeper::sweep1g( int group ) {


        return;
    }

    void MoCSweeper::initialize() {
        // There are better ways to do this, but for now, just start with 1.0
        flux_.fill(1.0);

        // Walk through the boundary conditions and initialize them the 1/4pi
        float_t val = 1.0/FPI;
        for( auto &group_rays: boundary_ ) {
            for( auto &plane_rays: group_rays ) {
                for( auto &angle_rays: plane_rays ) {
                    for( auto &ray: angle_rays ) {
                        ray.fw = val;
                        ray.bw = val;
                    }
                }
            }
        }
        return;
    }

    void MoCSweeper::calc_fission_source( float_t k, 
            MatrixX& fission_source ) const {

        float_t rkeff = 1.0/k;
        fission_source.fill(0.0);
        for( auto &xsr: xs_mesh_ ) {
            const auto& xsnf = xsr.xsmacnf();
            for(unsigned int ig=0; ig<xs_mesh_.n_grp(); ig++ ) {
                for( auto &ireg: xsr.reg() ) {
                    fission_source(ireg) += rkeff*xsnf[ig]*flux_(ig, ireg);
                }
            }
        }
        return;
    }

}
