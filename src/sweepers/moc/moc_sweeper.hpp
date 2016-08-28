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

#include <array>
#include "util/omp_guard.h"
#include "util/pugifwd.hpp"
#include "util/timers.hpp"
#include "core/angular_quadrature.hpp"
#include "core/boundary_condition.hpp"
#include "core/coarse_data.hpp"
#include "core/core_mesh.hpp"
#include "core/eigen_interface.hpp"
#include "core/exponential.hpp"
#include "core/transport_sweeper.hpp"
#include "core/xs_mesh.hpp"
#include "core/xs_mesh_homogenized.hpp"
#include "moc/ray_data.hpp"

namespace mocc {
namespace moc {
class MoCSweeper : public TransportSweeper {
public:
    MoCSweeper(const pugi::xml_node &input, const CoreMesh &mesh);

    ~MoCSweeper()
    {
    }

    virtual void sweep(int group) override;

    void initialize() override final;

    /**
     * \copydoc TransportSweeper::get_pin_flux_1g()
     *
     * The default \ref MeshTreatment for the \ref MoCSweeper and derived types
     * is \ref MeshTreatment::PLANE
     */
    void get_pin_flux_1g(
        int group, ArrayB1 &flux,
        MeshTreatment treatment = MeshTreatment::PLANE) const override final;

    /**
     * \copydoc TransportSweeper::set_pin_flux_1g()
     *
     * The default \ref MeshTreatment for \ref MoCSweeper is \ref
     * MeshTreatment::PLANE, which results in a pin-by-pin fine-mesh projection,
     * preserving the intra-pin flux shape for each pin. If passed \ref
     * MeshTreatment::PIN, an axial homogenization is performed first, and the
     * result is treated in the same way as MeshTreatment::PLANE.
     */
    real_t set_pin_flux_1g(
        int group, const ArrayB1 &pin_flux,
        MeshTreatment treatment = MeshTreatment::PLANE) override final;

    void output(H5Node &node) const override;

    void homogenize(CoarseData &data) const
    {
        throw EXCEPT("Not Implemented");
    }

    /**
     * \brief \copybrief TransportSweeper::update_incoming_flux()
     */
    void update_incoming_flux() override final;

    /**
     * \copybrief TransportSweeper::create_source()
     *
     * This mostly calls the base \ref TransportSweeper::create_source()
     * method, but also makes sure that the source is configured properly
     * for MoC.
     */
    virtual UP_Source_t
    create_source(const pugi::xml_node &input) const override final
    {
        auto source = TransportSweeper::create_source(input);
        return source;
    }

    /**
     * Return a copy of the sweeper's angular quadrature.
     */
    AngularQuadrature get_ang_quad() const
    {
        return ang_quad_;
    }

    SP_XSMeshHomogenized_t get_homogenized_xsmesh() override final
    {
        auto xsm = SP_XSMeshHomogenized_t(new XSMeshHomogenized(mesh_));
        xsm->set_flux(flux_);
        return xsm;
    }

    /**
     * \brief Apply a transverse leakage source
     *
     * This will apply the passed-in transverse leakage source to the
     * sweeper's source. If enabled and necessary, source splitting will be
     * used to enforce non-negativity on the external (non-self scatter)
     * source.
     */
    void apply_transverse_leakage(int group, const ArrayB1 &tl);

    void check_balance(int group) const;

protected:
    // Data
    Timer &timer_;
    Timer &timer_init_;
    Timer &timer_sweep_;
    const CoreMesh &mesh_;

    RayData rays_;

    // Multi-group, incoming boundary flux. One for each plane
    std::vector<BoundaryCondition> boundary_;
    // One-group, outgoing boundary flux
    std::vector<BoundaryCondition> boundary_out_;

    // Array of one group transport cross sections, including transverse
    // leakage splitting, if necessary
    ExpandedXS xstr_;

    // Reference to a one-group slice of flux_. This should be
    // default-constructed, so that it only references data in flux_
    ArrayB1 flux_1g_;

    // Subplane parameters. These come from the CoreMesh, ultimately through the
    // Assemblies. Each entry is the number of actual CoreMesh planes to bind
    // together into each macroplane
    VecI subplane_;

    // The upper bounds of each MoC plane. This is essentially the accumulation
    // of the entries in subplane_. This is useful for finding an MoC plane
    // index given a mesh z index using something like std::upper_bound().
    VecI subplane_bounds_;

    // Plane geometry IDs associated with each macroplane
    VecI macroplane_unique_ids_;

    // Vector containing the first region index in each macroplane
    VecI first_reg_macroplane_;

    // Number of FSRs in each MoC plane. These could be gleaned from the
    // CoreMesh, but storing them is just as easy
    VecI nreg_plane_;

    // The source splitting variable. This stores the degree by which to
    // alter the transport cross section for the current group
    ArrayB1 split_;

    // Number of inner iterations per group sweep
    unsigned int n_inner_;

    // Boundary condition enumeration
    std::array<Boundary, 6> bc_type_;

    // Exponential table
    Exponential_Linear<10000> exp_;

    bool dump_rays_;
    bool dump_fsr_flux_;
    bool gauss_seidel_boundary_;
    bool allow_splitting_;

    // Methods
    /**
     * \brief Return the MoC plane corresponding to the passed axial index
     */
    int moc_plane_index(int iz) const
    {
        return std::distance(subplane_bounds_.begin(),
                             std::upper_bound(subplane_bounds_.begin(),
                                              subplane_bounds_.end(), iz));
    }

#include "moc_sweeper_kernel.inc.hpp"

    template <class Function> void update_incoming_generic(Function f)
    {
        // There are probably more efficient ways to do this, but for now, just
        // loop over all of the rays, look up the appropriate surface from the
        // Mesh, and adjust the BC accordingly
        for (auto g : groups_) {
            int iplane = 0;
            for (auto plane_geom_id : mesh_.unique_plane_ids()) {
                auto &bc         = boundary_[iplane];
                const auto &rays = rays_[plane_geom_id];
                int iang         = 0;
                for (const auto &ang_rays : rays) {
                    int iang1 = iang;
                    int iang2 = ang_quad_.reverse(iang);

                    real_t *bc_fw = bc.get_boundary(g, iang1).second;
                    real_t *bc_bw = bc.get_boundary(g, iang2).second;

                    for (const auto &ray : ang_rays) {
                        int is1 = ray.cm_cell_fw();
                        int is2 = ray.cm_cell_bw();
                        int bc1 = ray.bc(0);
                        int bc2 = ray.bc(1);

                        bc_fw[bc1] = f(bc_fw[bc1], is1, g);
                        bc_bw[bc2] = f(bc_bw[bc2], is2, g);
                    } // rays

                    iang++;
                } // angles
                iplane++;
            } // planes
        }     // groups

        return;
    }
};
}
}
