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

#include "core/core_mesh.hpp"
#include "util/h5file.hpp"
#include "util/pugifwd.hpp"
#include "core/source.hpp"
#include "core/transport_sweeper.hpp"
#include "solver.hpp"

namespace mocc {
class FixedSourceSolver : public Solver {
public:
    /**
    * Initialize a FSS using an XML node and CoreMesh. This expects the passed
    * XML node to be a valid \<solver\> tag containing a relevant \<sweeper\>
    * tag, which is needed by the \ref TransportSweeperFactory() to generate a
    * \ref TransportSweeper.
    */
    FixedSourceSolver(const pugi::xml_node &input, const CoreMesh &mesh);

    ~FixedSourceSolver()
    {
    }

    /**
    * For now, there is no actual implementation of this method, since there is
    * no functionality for specifying a user-defined Source. In practice, the
    * FSS is driven via the \ref step() routine by the \ref EigenSolver.
    *
    * Ideally, this would solve a fixed source problem subject to the
    * configuration in the XML input. This can either be to some sort of
    * tolerance, or for a fixed number of group sweeps.
    */
    void solve();

    /**
    * Instructs the sweeper to store the old value of the flux, then performs a
    * group sweep.
    */
    void step();

    /**
    * Initialize the state of the FSS to start a new problem. For now this just
    * calls the same routine on the \ref TransportSweeper, which in turn
    * initializes the scalar flux, boundary conditions, etc. to some sort of
    * halfway-reasonable starting values.
    */
    void initialize()
    {
        sweeper_->initialize();
    }

    /**
     * Set the group-independent fission source. The group-dependent fission
     * source is calculated internally by the \ref Source object, typically at
     * the behest of an \ref EigenSolver
     */
    void set_fission_source(const ArrayB1 *fs)
    {
        assert((int)fs->size() == sweeper()->n_reg());
        fs_ = fs;
    }

    /**
     * Return the number of mesh regions.
     */
    unsigned int n_reg()
    {
        return sweeper_->n_reg();
    }

    /**
     * Return the number of energy groups
     */
    unsigned int n_group()
    {
        return ng_;
    }

    const TransportSweeper *sweeper() const
    {
        return sweeper_.get();
    }

    /**
     * Return a pointer to the the \ref TransportSweeper. Use with care.
     */
    TransportSweeper *sweeper()
    {
        return sweeper_.get();
    }

    void output(H5Node &node) const;

private:
    UP_Sweeper_t sweeper_;
    UP_Source_t source_;
    // Pointer to the group-independent fission source. Usually comes from an
    // eigenvalue solver, if present
    const ArrayB1 *fs_;
    size_t ng_;

    // Stuff that we should only need if we are doing a standalone FS solve
    bool fixed_source_;
    size_t max_iter_;
    real_t flux_tol_;
};
}
