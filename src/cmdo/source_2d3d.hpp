#pragma once

#include "error.hpp"
#include "source.hpp"

namespace mocc {
    /**
     * This class extends the Source class to provide an abstract representation
     * of two fused sources, which each treating one of the sweepers containted
     * in a PlaneSweeper_2D3D sweeper. Essentially, it is an opaque composition
     * of two sub-sources, which are both targeted by the various
     * sweeper-agnostic units in MOCC (e.g. EigenSolver and FixedSourceSolver)
     * by implementing fission() and self_scatter() methods which call the same
     * methods on the underlying Source objects. In effect, when EigenSolver
     * updates this source's fission source, both underlying sources get
     * updated. Likewise with the FixedSourceSolver updating the inscatter
     * source.
     *
     * The self_scatter() method, which is called by the sweeper itself, is not
     * implemented (or rather it will throw an error if it is called). This is
     * because the individual Sn and MoC sweepers should ultimately be assigned
     * their corresponding sub-sources, calling self_scatter() directly on those
     * instead.
     */
    class Source_2D3D: public Source {
    public:
        Source_2D3D( const MoCSweeper_2D3D &moc, const SnSweeper_CDD &sn ):
            Source( moc.n_reg(), &(moc.xs_mesh()), moc.flux() ),
            sn_source_( sn.n_reg(), &(sn.xs_mesh()), sn.flux() )
        {
            
        }

        void fission( const ArrayX &fs, int ig ) {
            Source::fission( fs, ig );
            sn_source_.fission( fs, ig );
        }

        void in_scatter( size_t ig ) {
            Source::in_scatter( ig );
            sn_source_.in_scatter( ig );
        }

        
        const Source* get_sn_source() const {
            return &sn_source_;
        }
    private:
        SnSource sn_source_;
    };
}
