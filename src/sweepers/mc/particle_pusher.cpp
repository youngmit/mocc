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

#include "particle.hpp"

namespace mocc {
ParticlePusher::ParticlePusher(const CoreMesh &mesh, const XSMesh &xs_mesh)
    : mesh_(mesh),
      xs_mesh_(xs_mesh),
      n_group_(xs_mesh.n_group()),
      do_implicit_capture_(true) {
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

void ParticlePusher::collide( Particle &p, int ixsreg ) {
    // Sample the type of interaction;
    const auto &xsreg = xs_mesh_[ixsreg];
    real_t r = RNG_MC.random();

    real_t xs_scat_out = xsreg.xsmacsc().out(p.group);
    real_t xs_cum = xs_scat_out;
    if(r < xs_cum) {
        // scatter. only isotropic for now
        // sample new energy
        real_t cdf_scale = 1.0/xs_scat_out;
        VecF scat_cdf;
        scat_cdf.reserve(n_group_);
        for( int iig=0; iig<n_group_; iig++ ) {
            scat_cdf.push_back(xsreg.xsmacsc().to(iig)[p.group]*cdf_scale);
        }
        assert(fp_equiv_ulp(1.0, scat_cdf.back()));

        p.group = RNG_MC.sample_cdf(scat_cdf);

        // sample new angle
        p.direction = Direction(RNG_MC.random(TWOPI), RNG_MC.random(-HPI, HPI));

        return;
    }

    xs_cum += xsreg.xsmacf(p.group);
    if(r < xs_cum) {
        // fission
        fission_bank_.push_back(p.location_global);
        return;
    }

    // collision
    if(!do_implicit_capture_) {
        p.alive = false;
    }

    return;
}

void ParticlePusher::simulate(Particle p) {
    // Figure out where the particle is
    auto location_info = mesh_.get_location_info(p.location_global,
                                                 p.direction);
    p.location = location_info.local_point;
    int ixsreg = xsmesh_regions_[location_info.ireg];

    real_t z_min = mesh_.z(location_info.pos.z);
    real_t z_max = mesh_.z(location_info.pos.z + 1);

    p.alive = true;

    while (p.alive) {
        // Determine distance to nearest surface
        auto d_to_surf = location_info.pm->distance_to_surface(p.location,
                                                               p.direction);

        // Determine distance to plane boundaries. If it is less than the
        // distance to surface, use it as the distance to surf and force a pin
        // intersection
        real_t d_to_plane =
            (p.direction.oz > 0.0)
                ? (z_max - p.location_global.z) / p.direction.oz
                : (p.location_global.z - z_min) / p.direction.oz;
        if (d_to_plane < d_to_surf.first) {
            d_to_surf.first = d_to_plane;
            d_to_surf.second =
                (p.direction.oz > 0.0) ? Surface::TOP : Surface::BOTTOM;
        }

        // Sample distance to collision
        real_t xstr = xs_mesh_[ixsreg].xsmactr(p.group);
        real_t d_to_collision = -std::log(1.0 - RNG_MC.random() * xstr) / xstr;

        if (d_to_collision < d_to_surf.first) {
            // Particle collided within the current region. Move particle to
            // collision site and handle interaction.
            this->collide( p, ixsreg );
        } else {
            // Particle reached a surface before colliding. Move to the
            // surface and re-sample distance to collision.
            if (d_to_surf.second == Surface::INTERNAL) {
                // Particle crossed an internal boundary in the pin.
                // Update its location and region index
                p.move(d_to_surf.first);
                location_info.ireg = location_info.pm->find_reg(p.location,
                                                                p.direction);
                ixsreg = xsmesh_regions_[location_info.ireg];
            } else {
                // Particle crossed a pin boundary. Move to neighboring pin,
                // handle boundary condition, etc.
                // Regardless of what happens, move the particle
                p.move(d_to_surf.first);

                // Figure out where we are again
                location_info = mesh_.get_location_info(p.location_global,
                                                        p.direction);
                z_min = mesh_.z(location_info.pos.z);
                z_max = mesh_.z(location_info.pos.z) + 1;
                p.location = location_info.local_point;
                // Check for domain boundary crossing
                if( location_info.surface == Surface::INTERNAL ) {
                    // We are still inside the domain.
                    // The update that we did above should be good enough that
                    // we don't need to do anything here
                } else {
                    // We are exiting a domain boundary. Handle the boundary
                    // condition.
                    auto bc = mesh_.boundary_condition(location_info.surface);
                    switch(bc) {
                    case Boundary::REFLECT:
                        p.direction.reflect(location_info.surface);
                        location_info = mesh_.get_location_info(
                                p.location_global, p.direction);
                        break;
                    case Boundary::VACUUM:
                        // Just kill the thing
                        p.alive = false;
                        break;
                    default:
                        throw EXCEPT("Unsupported boundary condition");
                    }
                }

                ixsreg = xsmesh_regions_[location_info.ireg];
            }
        }
    }
    return;
}

/**
 * This method operates by generating a full-blown Particle for each fission
 * site in the \ref FissionBank and calling \c this->simulate() for that
 * particle.
 */
void ParticlePusher::simulate(const FissionBank &bank) {
    // Clear the internal FissionBank to store new fission sites for this
    // cycle
    fission_bank_.clear();

    for (const auto &site : bank) {
        // Locate the particle in the problem geometry
        int ireg = mesh_.region_at_point(site);

        // Sample an energy group for the particle. If fissile, sample from
        // chi distribution, else sample uniformly
        int ixs = xsmesh_regions_[ireg];
        bool fissile = xs_mesh_[ixs].is_fissile();
        int ig = fissile ? RNG_MC.sample_cdf(xs_mesh_[ixs].chi_cdf())
                         : RNG_MC.random_int(n_group_);

        // Sample a direction
        real_t alpha = RNG_MC.random(TWOPI);
        real_t theta = RNG_MC.random(-1.0, 1.0);
        Direction d(alpha, theta);

        // Create the partile
        Particle p(site, d, ig);

        this->simulate(p);
    }
    return;
}

}  // namespace mocc
