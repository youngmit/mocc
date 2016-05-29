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

#include "particle_pusher.hpp"

#include <cmath>
#include <iostream>

#include "particle.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace mocc {
namespace mc {
ParticlePusher::ParticlePusher(const CoreMesh &mesh, const XSMesh &xs_mesh,
                               RNGSwarm &rng)
    : mesh_(mesh),
      xs_mesh_(xs_mesh),
      rng_(rng),
      n_group_(xs_mesh.n_group()),
      fission_bank_(mesh),
      do_implicit_capture_(false),
      scalar_flux_tally_(xs_mesh.n_group(), TallySpatial(mesh_.volumes()))
{
    // Build the map from mesh regions into the XS mesh
    xsmesh_regions_.reserve(mesh.n_reg());

    int ixs = 0;
    for (const auto &xsreg : xs_mesh_) {
        for (const auto &ireg : xsreg.reg()) {
            xsmesh_regions_[ireg] = ixs;
        }
        ixs++;
    }
    return;
}

void ParticlePusher::collide(Particle &p, int ixsreg)
{
    // Sample the type of interaction;
    const auto &xsreg = xs_mesh_[ixsreg];
    Reaction reaction =
        (Reaction)RNG_MC.sample_cdf(xsreg.reaction_cdf(p.group));

    switch (reaction) {
    case Reaction::SCATTER: {
        // scatter. only isotropic for now
        // sample new energy
        VecF cdf = xsreg.xsmacsc().out_cdf(p.group);
        p.group  = RNG_MC.sample_cdf(xsreg.xsmacsc().out_cdf(p.group));

        // sample new angle
        p.direction = Direction(RNG_MC.random(TWOPI), RNG_MC.random(-HPI, HPI));
    } break;

    case Reaction::FISSION:
        // fission
        // sample number of new particles to generate
        {
            real_t nu = xsreg.xsmacnf(p.group) / xsreg.xsmacf(p.group);
            assert(k_eff_ > 0.0);
            int n_fis = p.weight * nu / k_eff_ + RNG_MC.random();

            // Make new particles and push them onto the fission bank
            for (int i = 0; i < n_fis; i++) {
                int ig = RNG_MC.sample_cdf(xsreg.chi_cdf());
                Particle new_p(
                    p.location_global,
                    Direction(RNG_MC.random(TWOPI), RNG_MC.random(-HPI, HPI)),
                    ig);
                fission_bank_.push_back(new_p);
            }
        }
        p.alive = false;
        break;
    default:
        // capture
        p.alive = false;
    }

    return;
}

void ParticlePusher::simulate(Particle p)
{
    // Register this particle with the tallies
    k_tally_.add_weight(p.weight);
    for (auto &tally : scalar_flux_tally_) {
        tally.add_weight(p.weight);
    }

    p.alive = true;

    while (p.alive) {
        // Figure out where we are
        auto location_info =
            mesh_.get_location_info(p.location_global, p.direction);
        real_t z_min = mesh_.z(location_info.pos.z);
        real_t z_max = mesh_.z(location_info.pos.z) + 1;
        p.location   = location_info.local_point;
        int ireg     = location_info.reg_offset +
                   location_info.pm->find_reg(p.location, p.direction);
        assert(ireg >= 0);
        assert(ireg < mesh_.n_reg());
        int ixsreg = xsmesh_regions_[ireg];

        // Determine distance to nearest surface
        auto d_to_surf =
            location_info.pm->distance_to_surface(p.location, p.direction);
        // Determine distance to plane boundaries. If it is less than the
        // distance to surface, use it as the distance to surf and force a pin
        // intersection
        real_t d_to_plane =
            (p.direction.oz > 0.0)
                ? (z_max - p.location_global.z) / p.direction.oz
                : (z_min - p.location_global.z) / p.direction.oz;
        if (d_to_plane < d_to_surf.first) {
            d_to_surf.first = d_to_plane;
            d_to_surf.second =
                (p.direction.oz > 0.0) ? Surface::TOP : Surface::BOTTOM;
        }

        // Sample distance to collision
        real_t xstr           = xs_mesh_[ixsreg].xsmactr(p.group);
        real_t d_to_collision = -std::log(RNG_MC.random()) / xstr;

        real_t tl = std::min(d_to_collision, d_to_surf.first);

        // Contribute to track length-based tallies
        k_tally_.score(tl * p.weight * xs_mesh_[ixsreg].xsmacnf(p.group));

        scalar_flux_tally_[p.group].score(ireg, tl);

        if (d_to_collision < d_to_surf.first) {
            // Particle collided within the current region. Move particle to
            // collision site and handle interaction.
            p.move(d_to_collision);
            this->collide(p, ixsreg);
        }
        else {
            // Particle reached a surface before colliding. Move to the
            // surface and re-sample distance to collision.
            if (d_to_surf.second == Surface::INTERNAL) {
                // Particle crossed an internal boundary in the pin.
                // Update its location and region index
                p.move(d_to_surf.first);
                ireg = location_info.pm->find_reg(p.location, p.direction) +
                       location_info.reg_offset;
                ixsreg = xsmesh_regions_[ireg];
            }
            else {
                // Particle crossed a pin boundary. Move to neighboring pin,
                // handle boundary condition, etc.
                // Regardless of what happens, move the particle
                p.move(d_to_surf.first);

                // Check for domain boundary crossing
                Surface bound_surf =
                    mesh_.boundary_surface(p.location_global, p.direction);
                if (bound_surf != Surface::INTERNAL) {
                    // We are exiting a domain boundary. Handle the boundary
                    // condition.
                    auto bc = mesh_.boundary_condition(bound_surf);
                    switch (bc) {
                    case Boundary::REFLECT:
                        p.direction.reflect(bound_surf);
                        break;
                    case Boundary::VACUUM:
                        // Just kill the thing
                        p.alive = false;
                        break;
                    default:
                        throw EXCEPT("Unsupported boundary condition");
                    }
                }
            }
        } // collision or new region?
    }     // particle alive
    return;
}

/**
 * This method operates by generating a full-blown Particle for each fission
 * site in the \ref FissionBank and calling \c this->simulate() for that
 * particle.
 */
void ParticlePusher::simulate(const FissionBank &bank, real_t k_eff)
{
    // Clear the internal FissionBank to store new fission sites for this
    // cycle
    fission_bank_.clear();
    k_tally_.reset();

    k_eff_ = k_eff;

    for (const auto &p : bank) {
        this->simulate(p);
    }
    return;
}
} // namespace mc
} // namespace mocc
