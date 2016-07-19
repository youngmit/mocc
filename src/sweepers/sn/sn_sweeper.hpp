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

#include "util/pugifwd.hpp"
#include "util/timers.hpp"
#include "util/utils.hpp"
#include "core/angular_quadrature.hpp"
#include "core/boundary_condition.hpp"
#include "core/transport_sweeper.hpp"
#include "cell_worker.hpp"

namespace mocc {
namespace sn {
class SnSweeper : public TransportSweeper {
public:
    SnSweeper(const pugi::xml_node &input, const CoreMesh &mesh);

    void initialize()
    {
        flux_     = 1.0;
        flux_old_ = 1.0;
        bc_in_.initialize_scalar(1.0 / FPI);

        return;
    }

    void get_pin_flux_1g(int ig, ArrayB1 &flux) const
    {
        assert((int)flux.size() == n_reg_);

        flux = flux_(blitz::Range::all(), ig);

        return;
    }

    /**
     * \brief \copybrief TransportSweeper::update_incoming_flux()
     */
    void update_incoming_flux();

    ArrayB3 pin_powers() const;

    /**
     * Just copy the flux across, since no homogenization is necessary.
     */
    real_t set_pin_flux_1g(int group, const ArrayB1 &pin_flux)
    {
        assert((int)pin_flux.size() == n_reg_);

        real_t resid = 0.0;
        size_t i     = 0;
        for (auto &v : pin_flux) {
            real_t e = flux_(1, group) - v;
            resid += e * e;
            flux_((int)i, (int)group) = v;
            i++;
        }
        return std::sqrt(resid);
    }

    /**
     * \brief Re-assign the angular quadrature.
     */
    void set_ang_quad(AngularQuadrature ang_quad)
    {
        ang_quad_ = ang_quad;
        return;
    }

    SP_XSMeshHomogenized_t get_homogenized_xsmesh()
    {
        return std::static_pointer_cast<XSMeshHomogenized>(xs_mesh_);
    }

    ExpandedXS &expanded_xs() {
        return xstr_;
    }

    virtual void output(H5Node &node) const;

protected:
    Timer &timer_;
    Timer &timer_init_;
    Timer &timer_sweep_;
    Timer &timer_xsupdate_;
    const CoreMesh &mesh_;

    unsigned int n_inner_;

    // Boundary condition enumeration
    std::array<Boundary, 6> bc_type_;

    // One-group slice of flux_. Should be default-constructed, and assigned
    // slices using .reference()
    ArrayB1 flux_1g_;

    // Temporary storage of the current-group transport cross section
    ExpandedXS xstr_;

    // Incomming boundary condition
    BoundaryCondition bc_in_;

    // Outgoing boundary condition. Only difined for one group
    BoundaryCondition bc_out_;

    // Gauss-Seidel BC update?
    bool gs_boundary_;

    // Protected methods
    /**
     * \brief Grab data (XS, etc.) from one or more external files
     */
    void add_data(const pugi::xml_node &input);

    /**
     * \brief Check the neutron balance in all of the cells of the sweeper
     */
    void check_balance(int group) const;

    /**
     * This funtion template is used to permit flexibility in the incoming
     * flux update without incurring lots of code duplication, and withot
     * affecting performance. The loop over surfaces is somewhat involved,
     * and replicating it by hand for each flux update method would be
     * unmaintainable.
     *
     * The update for each surface is performed using a lambda function,
     * upon which this template is instantiated, which, under optimization,
     * should be inlined.
     */
    template <class Function> void update_incoming_generic(Function f)
    {
        for (auto g : groups_) {
            int iang = 0;
            for (auto ang : ang_quad_) {
                // X-normal
                {
                    real_t *face =
                        bc_in_.get_face(g, iang, Normal::X_NORM).second;
                    int ixx = (ang.ox > 0.0) ? 0 : mesh_.nx() - 1;
                    Surface upwind =
                        (ang.ox > 0.0) ? Surface::WEST : Surface::EAST;
                    int sense = (upwind == Surface::WEST) ? 1 : 0;
                    int i     = 0;
                    for (unsigned iz = 0; iz < mesh_.nz(); iz++) {
                        for (unsigned iy = 0; iy < mesh_.ny(); iy++) {
                            int icell =
                                mesh_.coarse_cell(Position(ixx, iy, iz));
                            int is = mesh_.coarse_surf(icell, upwind);

                            face[i] = f(face[i], is, g, sense);

                            i++;
                        }
                    }
                } // X-normal
                // Y-normal
                {
                    real_t *face =
                        bc_in_.get_face(g, iang, Normal::Y_NORM).second;
                    int iyy = (ang.oy > 0.0) ? 0 : mesh_.ny() - 1;
                    Surface upwind =
                        (ang.oy > 0.0) ? Surface::SOUTH : Surface::NORTH;
                    int sense = (upwind == Surface::SOUTH) ? 1 : 0;
                    int i     = 0;
                    for (unsigned iz = 0; iz < mesh_.nz(); iz++) {
                        for (unsigned ix = 0; ix < mesh_.nx(); ix++) {
                            int icell =
                                mesh_.coarse_cell(Position(ix, iyy, iz));
                            int is = mesh_.coarse_surf(icell, upwind);

                            face[i] = f(face[i], is, g, sense);

                            i++;
                        }
                    }
                } // Y-normal
                // Z-normal
                {
                    real_t *face =
                        bc_in_.get_face(g, iang, Normal::Z_NORM).second;
                    int izz = (ang.oz > 0.0) ? 0 : mesh_.nz() - 1;
                    Surface upwind =
                        (ang.oz > 0.0) ? Surface::BOTTOM : Surface::TOP;
                    int sense = (upwind == Surface::BOTTOM) ? 1 : 0;
                    int i     = 0;
                    for (unsigned iy = 0; iy < mesh_.ny(); iy++) {
                        for (unsigned ix = 0; ix < mesh_.nx(); ix++) {
                            int icell =
                                mesh_.coarse_cell(Position(ix, iy, izz));
                            int is = mesh_.coarse_surf(icell, upwind);

                            face[i] = f(face[i], is, g, sense);

                            i++;
                        }
                    }
                } // Z-normal

                iang++;
            } // angles
        }     // groups
    }

private:
};

typedef std::shared_ptr<SnSweeper> SP_SnSweeper_t;
typedef std::unique_ptr<SnSweeper> UP_SnSweeper_t;
}
}
