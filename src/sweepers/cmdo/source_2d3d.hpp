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

#pragma once

#include "util/error.hpp"
#include "core/source.hpp"

namespace mocc {
namespace cmdo {
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
class Source_2D3D : public SourceIsotropic {
public:
    Source_2D3D(const MoCSweeper_2D3D &moc, const TransportSweeper &sn)
        : SourceIsotropic(moc.n_reg(), &(moc.xs_mesh()), moc.flux()),
          mesh_(moc.mesh()),
          sn_source_(sn.n_reg(), &(sn.xs_mesh()), sn.flux())
    {
        return;
    }

    /**
     * Replace the standard group initializer with a call to the base type
     * and the Sn source inside.
     */
    void initialize_group(int group)
    {
        Source::initialize_group(group);
        sn_source_.initialize_group(group);
        return;
    }

    /**
     * Replaces the standard fission source calculation with a delegation to the
     * base \ref Source::fission() routine for MoC, a possible homogenization of
     * the FM fission source to the Sn mesh, and a call to \ref
     * Source::fission() on the Sn source object with the homogenized fission
     * source.
     *
     * There are two possibilities for how the fission source might be defined:
     * either on the fine mesh or coarse, Sn mesh. For now it should be
     * sufficient to use the size of the \p fs array to decide which we are
     * getting.
     */
    void fission(const ArrayB1 &fs, int ig)
    {
        assert((int)fs.size() == (n_reg_ + sn_source_.n_reg()));

        // Alias the appropriate regions of the passed fission source
        ArrayB1 sn_fission_source(fs(blitz::Range(0, sn_source_.n_reg() - 1)));
        ArrayB1 moc_fission_source(
            fs(blitz::Range(sn_source_.n_reg(), blitz::toEnd)));

        sn_source_.fission(sn_fission_source, ig);
        Source::fission(moc_fission_source, ig);
    }

    void in_scatter(size_t ig)
    {
        Source::in_scatter(ig);
        sn_source_.in_scatter(ig);
    }

    Source *get_sn_source()
    {
        return &sn_source_;
    }

private:
    const CoreMesh &mesh_;
    SourceIsotropic sn_source_;
};
} // namespace cmdo
} // namespace mocc
