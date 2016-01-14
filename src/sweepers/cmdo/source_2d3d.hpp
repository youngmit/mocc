#pragma once

#include "core/error.hpp"
#include "core/source.hpp"

namespace mocc {
    /**
     * This class extends the Source class to provide an abstract representation
     * of two fused sources, with each treating one of the sweepers containted
     * in a PlaneSweeper_2D3D sweeper. Essentially, it is an opaque composition
     * of two sub-sources, which are both targeted by the various
     * sweeper-agnostic units in MOCC (e.g. EigenSolver and FixedSourceSolver)
     * by implementing fission() and self_scatter() methods which call the same
     * methods on the underlying Source objects, and perform whatever
     * homogenization operations are needed. In effect, when EigenSolver updates
     * this source's fission source, both underlying sources get updated.
     * Likewise with the FixedSourceSolver updating the inscatter source.
     *
     * It should also be noted that the base Source class from which this
     * inherits is used as the MoC source, and all of the base data members
     * should be treated as though they apply to the MoC sweeper.
     *
     * The self_scatter() method, which is called by the sweeper itself, is not
     * implemented (or rather it will throw an error if it is called). This is
     * because the individual Sn and MoC sweepers should ultimately be assigned
     * their corresponding sub-sources, calling self_scatter() directly on those
     * instead.
     */
    class Source_2D3D: public SourceIsotropic {
    public:
        Source_2D3D( const MoCSweeper_2D3D &moc, const TransportSweeper &sn ):
            SourceIsotropic( moc.n_reg(), &(moc.xs_mesh()), moc.flux() ),
            mesh_(moc.mesh()),
            sn_source_( sn.n_reg(), &(sn.xs_mesh()), sn.flux() )
        {
            this->set_scale_transport(true);
        }

        /**
         * Replace the standard group initializer with a call to the base type
         * and the Sn source inside.
         */
        void initialize_group( int group ) {
            Source::initialize_group( group );
            sn_source_.initialize_group( group );
            return;
        }

        /**
         * Replaces the standard fission source calculation with a delegation to
         * the base \ref Source::fission() routine for MoC, a possible
         * homogenization of the FM fission source to the Sn mesh, and a call to
         * \ref Source::fission() on the Sn source object with the homogenized
         * fissions source.
         *
         * There are two possibilities for how the fission source might be
         * defined: either on the fine mesh or coarse, Sn mesh. For now it
         * should be sufficient to use the size of the \p fs array to decide
         * which we are getting.
         */
        void fission( const ArrayF &fs, int ig ) {
            assert( (fs.size() == n_reg_) ||
                    (fs.size() == sn_source_.n_reg()) );

            Source::fission( fs, ig );

            if( fs.size() == n_reg_ ) {
                // We need to homogenize the fission source to the Sn mesh
                ArrayF sn_fs(sn_source_.n_reg());
                int ireg_fsr = 0;
                int ipin = 0;
                for( const auto &pin: mesh_ ) {
                    auto pos = mesh_.pin_position(ipin);
                    int ireg = mesh_.index_lex( pos );
                    for( const auto v: pin->vols() ) {
                        sn_fs[ireg] += v*fs[ireg_fsr];
                        ireg_fsr++;
                    }
                    sn_fs[ireg] /= pin->vol();
                    ipin++;
                }

                sn_source_.fission( sn_fs, ig );
            } else {
                // Let's assume the fissions source is defined on the Sn mesh,
                // and set it directly
                sn_source_.fission( fs, ig );
            }
        }

        void in_scatter( size_t ig ) {
            Source::in_scatter( ig );
            sn_source_.in_scatter( ig );
        }


        Source* get_sn_source() {
            return &sn_source_;
        }
    private:
        const CoreMesh& mesh_;
        SourceIsotropic sn_source_;
    };
}
