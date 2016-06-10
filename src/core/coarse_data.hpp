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

#include <memory>
#include <vector>

#include <blitz/array.h>

#include "core/blitz_typedefs.hpp"
#include "core/eigen_interface.hpp"
#include "core/global_config.hpp"
#include "core/mesh.hpp"

namespace mocc {
/**
 * CoarseData stores the data needed to do CMFD. Coarse surface currents,
 * fluxes, etc.
 *
 * \todo the storage order is reversed from what it should be. When things
 * settle down some, do some profiling and swap
 */
struct CoarseData {
public:
    CoarseData(const Mesh &mesh, size_t ngroup);

    /**
     * \brief Signal to other clients of the \ref CoarseData that data have
     * been explicitly defined for the radial-facing surfaces (X- and
     * Y-normal).
     *
     * Upon construction, this is set to \c false. When a sweeper or other
     * client sets values on the \ref CoarseData, it should set it to true
     * so that other clients know that they can use it, via the \ref
     * has_radial_data() method.
     *
     * For example, the \ref CMFD solver should not try to calculate D-hats
     * for surfaces unless currents have been supplied by a transport
     * sweeper. In the context of a 2-D MoC sweeper, which never sets axial
     * currents, the \ref CMFD solver should never calculate D-hats, even
     * though the axial currents may be non-zero (since the \ref CMFD solver
     * updats the currents after a solve).
     */
    void set_has_radial_data(bool has)
    {
        has_data_radial_ = has;
        return;
    }

    /**
     * \brief Signal to other clients of the \ref CoarseData that axial
     * (Z-normal) currents have been explicitly defined.
     *
     * This is similar to \ref set_has_radial_data(), but for the axial
     * surfaces. This should be set \c true after a 3-D sweeper has
     * calculated currents.
     */
    void set_has_axial_data(bool has)
    {
        has_data_axial_ = has;
        return;
    }

    /**
     * \brief Return whether data have been explicitly defined for axial
     * surfaces.
     */
    bool has_axial_data() const
    {
        return has_data_axial_;
    }

    /**
     * \brief Return whether data have been explicitly defined for radial
     * surfaces.
     */
    bool has_radial_data() const
    {
        return has_data_radial_;
    }

    /**
     * \brief Return whether previous-iteration values for partial currents
     * are available.
     *
     * This starts as \c false at construction time, and is set to \c true
     * immediately after the first \ref CMFD solve, or similar. This is
     * necessary, since the logic for tasks like updating incomming flux is
     * different if old values are available. See the various
     * implementations of \ref TransportSweeper::update_incoming_flux().
     */
    bool has_old_partial() const
    {
        return has_old_partial_;
    }

    /**
     * \brief Signal to other clients of the \ref CoarseData that it has
     * previous-iteration values for partial current.
     */
    void set_has_old_partial(bool has)
    {
        has_old_partial_ = has;
        return;
    }

    /**
     * \brief Zero out all of the data associated with the given group.
     *
     * This is typically used immediately before invoking a sweep procedure
     * that will calculate new data.
     *
     * \note This will zero out all of the data, including radial and axial
     * surfaces. It is therefore best suited for use with 3-D sweepers. Most
     * 2-D sweepers will want to use the 2-D version, \ref
     * zero_data_radial().
     */
    void zero_data(int group, bool zero_partial = false);

    /**
     * \brief Zero out the data on the radial-normal surfaces for a given
     * group.
     *
     * This is the 2-D version of \ref zero_data(). It zeros out the X- and
     * Y-normal surfaces, but leaves data for the other surfaces untouched.
     */
    void zero_data_radial(int group, bool zero_partial = false);

    ArrayB2 current;
    ArrayB2 surface_flux;
    blitz::Array<std::array<real_t, 2>, 2> partial_current;
    blitz::Array<std::array<real_t, 2>, 2> partial_current_old;
    ArrayB2 flux;
    ArrayB2 old_flux;

private:
    int n_group_;
    const Mesh &mesh_;
    bool has_data_radial_;
    bool has_data_axial_;
    bool has_old_partial_;
};

typedef std::shared_ptr<CoarseData> SP_CoarseData_t;
}
