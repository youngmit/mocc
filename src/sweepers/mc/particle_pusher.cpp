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
#include <omp.h>

#include "core/blitz_typedefs.hpp"
#include "core/utils.hpp"

#include "particle.hpp"

using std::cout;
using std::endl;
using std::cin;

namespace {
using namespace mocc;
const real_t BUMP = 1.0e-11;
}

namespace mocc {
namespace mc {

// Extern hack needed to get GCC to allow threadprivate instance of non-POD type
extern RNG_LCG RNG;
#pragma omp threadprivate(RNG)
RNG_LCG RNG;

ParticlePusher::ParticlePusher(const CoreMesh &mesh, const XSMesh &xs_mesh)
    : mesh_(mesh),
      xs_mesh_(xs_mesh),
      n_group_(xs_mesh.n_group()),
      fission_bank_(mesh),
      do_implicit_capture_(false),
      seed_(1),
      scalar_flux_tally_(xs_mesh.n_group(), TallySpatial(mesh_.volumes())),
      id_offset_(0),
      n_cycles_(0),
      print_particles_(false)
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
    Reaction reaction = (Reaction)RNG.sample_cdf(xsreg.reaction_cdf(p.group));
    real_t k_score    = p.weight * xs_mesh_[ixsreg].xsmacnf(p.group) /
                     xs_mesh_[ixsreg].xsmactr(p.group);
    k_tally_collision_.score(k_score);

    switch (reaction) {
    case Reaction::SCATTER: {
        if (print_particles_)
            cout << "scatter from group " << p.group << endl;
        // scatter. only isotropic for now
        // sample new energy
        p.group = RNG.sample_cdf(xsreg.xsmacsc().out_cdf(p.group));
        if (print_particles_) {
            cout << "New group: " << p.group << endl;
        }

        // sample new angle
        p.direction = Direction(RNG.random(TWOPI), RNG.random(PI));
        if (print_particles_) {
            cout << "New angle: " << p.direction << endl;
        }
    } break;

    case Reaction::FISSION:
        if (print_particles_) {
            cout << "fission at " << p.location_global.x << " "
                 << p.location_global.y << " " << p.location_global.z << endl;
        }
        // fission
        // sample number of new particles to generate
        {
            real_t nu = xsreg.xsmacnf(p.group) / xsreg.xsmacf(p.group);
            assert(k_eff_ > 0.0);
            int n_fis = p.weight * nu / k_eff_ + RNG.random();

            // Make new particles and push them onto the fission bank
            for (int i = 0; i < n_fis; i++) {
                int ig = RNG.sample_cdf(xsreg.chi_cdf());
                Particle new_p(p.location_global,
                               Direction(RNG.random(TWOPI), RNG.random(PI)), ig,
                               p.id);
                fission_bank_.push_back(new_p);
            }
        }
        p.alive = false;
        break;
    default:
        if (print_particles_) {
            cout << "capture" << endl;
        }
        // capture
        p.alive = false;
    }

    return;
}

void ParticlePusher::simulate(Particle p)
{
    if (print_particles_)
        cout << endl << "NEW PARTICLE:" << endl;
    // Register this particle with the tallies
    k_tally_.add_weight(p.weight);
    k_tally_collision_.add_weight(p.weight);
    for (auto &tally : scalar_flux_tally_) {
        tally.add_weight(p.weight);
    }
    // Figure out where we are
    auto location_info =
        mesh_.get_location_info(p.location_global, p.direction);
    real_t z_min = mesh_.z(location_info.pos.z);
    real_t z_max = mesh_.z(location_info.pos.z + 1);
    p.location   = location_info.local_point;
    int ireg     = location_info.reg_offset +
               location_info.pm->find_reg(p.location, p.direction);
    assert(ireg >= 0);
    assert(ireg < (int)mesh_.n_reg());
    int ixsreg = xsmesh_regions_[ireg];

    p.alive = true;

    while (p.alive) {
        if (print_particles_) {
            cout << "Where we are now:" << endl;
            cout << p << endl;
            cout << "ireg/xsreg: " << ireg << " " << ixsreg << endl;
        }
        // Determine distance to nearest surface
        auto d_to_surf =
            location_info.pm->distance_to_surface(p.location, p.direction);
        // Determine distance to plane boundaries. If it is less than the
        // distance to surface, use it as the distance to surf and force a pin
        // intersection
        real_t d_to_plane =
            (p.direction.oz > 0.0)
                ? std::max(0.0, (z_max - p.location_global.z) / p.direction.oz)
                : std::max(0.0, (z_min - p.location_global.z) / p.direction.oz);
        if (d_to_plane < d_to_surf.first) {
            d_to_surf.first = d_to_plane;
            d_to_surf.second =
                (p.direction.oz > 0.0) ? Surface::TOP : Surface::BOTTOM;
        }

        // Sample distance to collision
        real_t xstr           = xs_mesh_[ixsreg].xsmactr(p.group);
        real_t d_to_collision = -std::log(RNG.random()) / xstr;

        if (print_particles_) {
            cout << "distance to surface/collision: " << d_to_surf.first << " "
                 << d_to_collision << endl;
        }

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
                p.move(d_to_surf.first + BUMP);
                ireg = location_info.pm->find_reg(p.location, p.direction) +
                       location_info.reg_offset;
                ixsreg = xsmesh_regions_[ireg];
            }
            else {
                // Particle crossed a pin boundary. Move to neighboring pin,
                // handle boundary condition, etc.
                // Regardless of what happens, move the particle
                p.move(d_to_surf.first + BUMP);

                if (print_particles_) {
                    cout << "particle after move to surf:" << endl;
                    cout << p << endl;
                }

                // Check for domain boundary crossing
                auto bound_surf =
                    mesh_.boundary_surface(p.location_global, p.direction);
                bool reflected = false;
                for (const auto &b : bound_surf) {
                    if (print_particles_) {
                        cout << b << endl;
                    }
                    if ((b != Surface::INTERNAL) && (p.alive)) {
                        // We are exiting a domain boundary. Handle the boundary
                        // condition.
                        auto bc = mesh_.boundary_condition(b);
                        switch (bc) {
                        case Boundary::REFLECT:
                            // Move the particle back into the domain so it's
                            // not
                            // floating in limbo
                            reflected = true;
                            p.direction.reflect(b);
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

                if (reflected) {
                    p.move(2.0 * BUMP);
                    if (print_particles_) {
                        cout << "Particle after reflection and move back:"
                             << endl;
                        cout << p << endl;
                    }
                }

                // If the particle is still alive, relocate it
                if (p.alive) {
                    location_info =
                        mesh_.get_location_info(p.location_global, p.direction);
                    z_min      = mesh_.z(location_info.pos.z);
                    z_max      = mesh_.z(location_info.pos.z + 1);
                    p.location = location_info.local_point;
                    ireg       = location_info.reg_offset +
                           location_info.pm->find_reg(p.location, p.direction);
                    assert(ireg >= 0);
                    assert(ireg < (int)mesh_.n_reg());
                    ixsreg = xsmesh_regions_[ireg];
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
    k_tally_collision_.reset();

    k_eff_ = k_eff;

    print_particles_ = false;

#pragma omp parallel
    {
        unsigned np = bank.size();
#pragma omp for
        for (unsigned ip = 0; ip < np; ip++) {
            RNG.set_seed(seed_);
            RNG.jump_ahead((bank[ip].id + id_offset_) * 10000);
            this->simulate(bank[ip]);
        }

        n_cycles_++;
    } // OMP Parallel

    id_offset_ += bank.size();
    return;
}

void ParticlePusher::output(H5Node &node) const
{
    auto dims = mesh_.dimensions();
    std::reverse(dims.begin(), dims.end());

    node.write("ng", n_group_);
    node.write("eubounds", xs_mesh_.eubounds(), VecI(1, n_group_));

    ArrayB2 flux_mg(n_group_, mesh_.n_pin());
    ArrayB2 stdev_mg(n_group_, mesh_.n_pin());

    node.create_group("flux");
    for (int ig = 0; ig < (int)scalar_flux_tally_.size(); ig++) {
        auto flux_result = scalar_flux_tally_[ig].get_homogenized(mesh_);
        int ipin         = 0;
        for (const auto &v : flux_result) {
            flux_mg(ig, ipin)  = v.first;
            stdev_mg(ig, ipin) = v.second;
            ipin++;
        }
    }

    real_t f = Normalize(flux_mg.begin(), flux_mg.end());
    Scale(stdev_mg.begin(), stdev_mg.end(), f);

    for (int ig = 0; ig < (int)scalar_flux_tally_.size(); ig++) {
        std::stringstream path;
        path << "flux/" << std::setfill('0') << std::setw(3) << ig + 1;

        node.write(path.str(), flux_mg(ig, blitz::Range::all()), dims);
        path << "_stdev";
        node.write(path.str(), stdev_mg(ig, blitz::Range::all()), dims);
    }
    return;
}
} // namespace mc
} // namespace mocc
