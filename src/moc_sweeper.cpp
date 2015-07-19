#include "moc_sweeper.hpp"

#include "error.hpp"

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

        return;
    }

    void MoCSweeper::sweep( int group ) {
        return;
    }

    void MoCSweeper::sweep1g( int group ) {


        return;
    }

    void MoCSweeper::initialize() {
        // There are better ways to do this, but for now, just start with 1.0
        flux_.fill(1.0);
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
