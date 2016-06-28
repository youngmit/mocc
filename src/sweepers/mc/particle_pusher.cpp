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
#include "util/blitz_typedefs.hpp"
#include "util/omp_guard.h"
#include "util/utils.hpp"
#include "particle.hpp"

namespace {
using namespace mocc;
// Given an array of two points and a Particle, determine the distance from the
// particle to the box formed by the
real_t inline distance_to_pin(std::array<Point3, 2> bounds,
                              const mc::Particle &p)
{
    real_t dist  = std::numeric_limits<real_t>::max();
    real_t other = 0.0;

    other = (p.direction.ox > 0.0)
                ? (bounds[1].x - p.location_global.x) / p.direction.ox
                : (bounds[0].x - p.location_global.x) / p.direction.ox;
    other = other > 0.0 ? other : std::numeric_limits<real_t>::max();
    dist  = std::min(other, dist);

    other = (p.direction.oy > 0.0)
                ? (bounds[1].y - p.location_global.y) / p.direction.oy
                : (bounds[0].y - p.location_global.y) / p.direction.oy;
    other = other > 0.0 ? other : std::numeric_limits<real_t>::max();
    dist  = std::min(other, dist);

    other = (p.direction.oz > 0.0)
                ? (bounds[1].z - p.location_global.z) / p.direction.oz
                : (bounds[0].z - p.location_global.z) / p.direction.oz;
    other = other > 0.0 ? other : std::numeric_limits<real_t>::max();
    dist  = std::min(other, dist);

    assert(dist < std::numeric_limits<real_t>::max());
    assert(dist >= 0.0);
    return dist;
}
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
      volumes_(mesh.volumes(MeshTreatment::TRUE)),
      n_group_(xs_mesh.n_group()),
      fission_bank_(mesh),
      do_implicit_capture_(false),
      seed_(1),
      scalar_flux_tally_(xs_mesh.n_group(),
                         TallySpatial(mesh_.coarse_volume())),
      fine_flux_tally_(xs_mesh.n_group(), TallySpatial(volumes_)),
      fine_flux_col_tally_(xs_mesh.n_group(), TallySpatial(volumes_)),
      pin_power_tally_(mesh_.coarse_volume()),
      id_offset_(0),
      n_cycles_(0),
      print_particles_(false)
{
    // Build the map from mesh regions into the XS mesh
    xsmesh_regions_.resize(mesh.n_reg(MeshTreatment::TRUE), -1);

    int ixs = 0;
    for (const auto &xsreg : xs_mesh_) {
        for (const auto &ireg : xsreg.reg()) {
            xsmesh_regions_[ireg] = ixs;
        }
        ixs++;
    }
    return;
}

void ParticlePusher::collide(Particle &p)
{
    bool print = print_particles_;
    // print      = true;
    // Sample the type of interaction;
    const auto &xsreg = xs_mesh_[p.ixsreg];
    if (print) {
        std::cout << "COLLISION" << std::endl;
        std::cout << p << std::endl;
        std::cout << "xsregion: " << p.ixsreg << std::endl;
        std::cout << "reaction chance: ";
        for (const auto &v : xsreg.reaction_cdf(p.group)) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
    Reaction reaction = (Reaction)RNG.sample_cdf(xsreg.reaction_cdf(p.group));
    real_t k_score = p.weight * xsreg.xsmacnf(p.group) / xsreg.xsmactr(p.group);
    k_tally_col_.score(k_score);
    fine_flux_col_tally_[p.group].score(p.ireg,
                                        p.weight / xsreg.xsmactr(p.group));

    if (reaction == Reaction::SCATTER) {
        if (print) {
            std::cout << "scatter from group " << p.group << std::endl;
        }
        // scatter. only isotropic for now
        // sample new energy
        p.group = RNG.sample_cdf(xsreg.xsmacsc().out_cdf(p.group));
        if (print) {
            std::cout << "New group: " << p.group << std::endl;
        }

        // sample new angle
        p.direction = Direction::Isotropic(RNG.random(), RNG.random());
        if (print) {
            std::cout << "New angle: " << p.direction << std::endl;
        }
    } else if (reaction == Reaction::FISSION) {
        if (print) {
            std::cout << "fission at " << p.location_global.x << " "
                 << p.location_global.y << " " << p.location_global.z << std::endl;
        }
        // fission
        // sample number of new particles to generate
        real_t nu = xsreg.xsmacnf(p.group) / xsreg.xsmacf(p.group);
        int n_fis = (p.weight * nu) + RNG.random();

        // Make new particles and push them onto the fission bank
        for (int i = 0; i < n_fis; i++) {
            int ig = RNG.sample_cdf(xsreg.chi_cdf());
            Particle new_p(p.location_global,
                           Direction::Isotropic(RNG.random(), RNG.random()), ig,
                           p.id);
            fission_bank_.push_back(new_p);
        }
        p.alive = false;
    } else {
        if (print) {
            std::cout << "capture" << std::endl;
        }
        // capture
        p.alive = false;
    }

    return;
}

void ParticlePusher::simulate(Particle p, bool tally)
{
    bool print = print_particles_;

    RNG.set_seed(seed_);
    RNG.jump_ahead((p.id + id_offset_) * 10000);

    // Register this particle with the tallies
    k_tally_tl_.add_weight(p.weight);
    k_tally_col_.add_weight(p.weight);
    for (auto &tally : scalar_flux_tally_) {
        tally.add_weight(p.weight);
    }
    for (auto &tally : fine_flux_tally_) {
        tally.add_weight(p.weight);
    }
    for (auto &tally : fine_flux_col_tally_) {
        tally.add_weight(p.weight);
    }
    pin_power_tally_.add_weight(p.weight);

    // Figure out where we are
    auto location_info =
        mesh_.get_location_info(p.location_global, p.direction);
    p.location = location_info.local_point;
    p.ireg     = location_info.reg_offset +
             location_info.pm->find_reg(p.location, p.direction);
    p.pin_position  = location_info.pos;
    int ipin_coarse = mesh_.coarse_cell(location_info.pos);
    if (print) {
        std::cout << std::endl << "NEW PARTICLE:" << std::endl;
        std::cout << p << std::endl;
    }
    assert(ipin_coarse >= 0);
    assert(ipin_coarse < (int)mesh_.n_pin());
    assert(p.ireg >= 0);
    assert(p.ireg < (int)mesh_.n_reg(MeshTreatment::TRUE));
    p.ixsreg = xsmesh_regions_[p.ireg];
    assert(p.ixsreg >= 0);
    assert(p.ixsreg < (int)xs_mesh_.size());

    p.alive = true;

    while (p.alive) {
        const XSMeshRegion &xsreg = xs_mesh_[p.ixsreg];
        real_t xstr               = xsreg.xsmactr(p.group);
        real_t d_to_collision     = -std::log(RNG.random()) / xstr;

        // Determine distance to nearest surface
        auto d_to_surf = location_info.pm->distance_to_surface(
            p.location, p.direction, p.coincident);
        if (print) {
            std::cout << "Where we are now:" << std::endl;
            std::cout << p << std::endl;
            std::cout << "ireg/xsreg: " << p.ireg << " " << p.ixsreg << std::endl;
            std::cout << "Distance to internal pin surf: " << d_to_surf.first << " "
                 << d_to_surf.second << std::endl;
        }
        // Determine distance to plane boundaries. If it is less than the
        // distance to surface, use it as the distance to surf and force a pin
        // intersection
        real_t d_to_pin = distance_to_pin(location_info.pin_boundary, p);
        if (d_to_pin < d_to_surf.first) {
            d_to_surf.first  = d_to_pin;
            d_to_surf.second = true;
        }

        if (print) {
            std::cout << "distance to surface/collision: " << d_to_surf.first << " "
                 << d_to_surf.second << " " << d_to_collision << std::endl;
        }

        real_t tl = std::min(d_to_collision, d_to_surf.first);

        // Contribute to track length-based tallies
        k_tally_tl_.score(tl * p.weight * xsreg.xsmacnf(p.group));
        pin_power_tally_.score(ipin_coarse,
                               tl * p.weight * xsreg.xsmacf(p.group));
        scalar_flux_tally_[p.group].score(ipin_coarse, tl * p.weight);
        fine_flux_tally_[p.group].score(p.ireg, tl * p.weight);

        if (d_to_collision < d_to_surf.first) {
            // Particle collided within the current region. Move particle to
            // collision site and handle interaction.
            p.move(d_to_collision);
            p.coincident = -1;
            if (print) {
                std::cout << "particle at collision site:" << std::endl;
                std::cout << p << std::endl;
            }
            this->collide(p);
        } else {
            // Particle reached a surface before colliding. Move to the
            // surface and re-sample distance to collision.
            if (!d_to_surf.second) {
                // Particle crossed an internal boundary in the pin.
                // Update its location and region index
                p.move(d_to_surf.first);
                p.ireg = location_info.pm->find_reg(p.location, p.direction) +
                         location_info.reg_offset;
                assert(p.ireg >= 0);
                assert(p.ireg < (int)mesh_.n_reg(MeshTreatment::TRUE));
                p.ixsreg = xsmesh_regions_[p.ireg];

                assert(p.ixsreg >= 0);
                assert(p.ixsreg < (int)xs_mesh_.size());
            } else {
                // Particle crossed a pin boundary. Move to neighboring pin,
                // handle boundary condition, etc.
                // Regardless of what happens, move the particle
                p.coincident = -1;
                p.move(d_to_surf.first);

                if (print) {
                    std::cout << "particle after move to surf:" << std::endl;
                    std::cout << p << std::endl;
                }

                // Check for domain boundary crossing
                auto bound_surf =
                    mesh_.boundary_surface(p.location_global, p.direction);
                bool reflected = false;
                for (const auto &b : bound_surf) {
                    if (print) {
                        std::cout << b << std::endl;
                    }
                    if ((b != Surface::INTERNAL) && (p.alive)) {
                        // We are exiting a domain boundary. Handle the boundary
                        // condition.
                        auto bc = mesh_.boundary_condition(b);
                        switch (bc) {
                        case Boundary::REFLECT:
                            // Move the particle back into the domain so it's
                            // not floating in limbo
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
                    p.move(BUMP);
                    if (print) {
                        std::cout << "Particle after reflection and move back:"
                             << std::endl;
                        std::cout << p << std::endl;
                    }
                }

                // If the particle is still alive, relocate it
                if (p.alive) {
                    location_info =
                        mesh_.get_location_info(p.location_global, p.direction);
                    p.location = location_info.local_point;
                    p.ireg =
                        location_info.reg_offset +
                        location_info.pm->find_reg(p.location, p.direction);
                    p.pin_position = location_info.pos;
                    ipin_coarse    = mesh_.coarse_cell(location_info.pos);
                    assert(ipin_coarse >= 0);
                    assert(ipin_coarse < (int)mesh_.n_pin());
                    assert(p.ireg >= 0);
                    assert(p.ireg < (int)mesh_.n_reg(MeshTreatment::TRUE));
                    p.ixsreg = xsmesh_regions_[p.ireg];
                    if (p.ixsreg >= (int)xs_mesh_.size()) {
                        std::cout << p;
                        std::cout << std::endl;
                    }
                    assert(p.ixsreg >= 0);
                    assert(p.ixsreg < (int)xs_mesh_.size());
                }
            }
        } // collision or new region?
    }     // particle alive

    if (tally) {
        this->commit_tallies();
    }

    return;
}

/**
 * \brief Simulate all particles in a \ref FissionBank, stashing statistics at
 * the end.
 */
void ParticlePusher::simulate(const FissionBank &bank, real_t k_eff)
{
    // Clear the internal FissionBank to store new fission sites for this
    // cycle
    fission_bank_.clear();

    // not even used for now. Just letting the fission bank grow, then resizing
    // at the end
    k_eff_ = k_eff;

    print_particles_ = false;

#pragma omp parallel
    {
        unsigned np = bank.size();
#pragma omp for
        for (unsigned ip = 0; ip < np; ip++) {
            this->simulate(bank[ip]);
        }

    } // OMP Parallel

    k_tally_analog_.score((double)fission_bank_.size() / bank.size());
    k_tally_analog_.add_weight(1.0);

    this->commit_tallies();

    n_cycles_++;

    id_offset_ += bank.size();
    return;
}

void ParticlePusher::output(H5Node &node) const
{
    auto dims = mesh_.dimensions();
    std::reverse(dims.begin(), dims.end());

    node.write("ng", n_group_);
    node.write("eubounds", xs_mesh_.eubounds(), VecI(1, n_group_));

    // Coarse flux tallies
    {
        ArrayB2 flux_mg(xs_mesh_.n_group(), mesh_.n_pin());
        ArrayB2 stdev_mg(xs_mesh_.n_group(), mesh_.n_pin());

        auto g = node.create_group("flux");
        for (int ig = 0; ig < (int)scalar_flux_tally_.size(); ig++) {
            auto flux_result = scalar_flux_tally_[ig].get();
            int ipin         = 0;
            for (const auto &v : flux_result) {
                flux_mg(ig, ipin)  = v.first;
                stdev_mg(ig, ipin) = v.second;
                ipin++;
            }
        }

        /// \todo make flux normalization optional
        // Normalize(flux_mg.begin(), flux_mg.end());

        for (int ig = 0; ig < (int)scalar_flux_tally_.size(); ig++) {
            std::stringstream path;
            path << std::setfill('0') << std::setw(3) << ig + 1;

            g.write(path.str(), flux_mg(ig, blitz::Range::all()), dims);
            path << "_stdev";
            g.write(path.str(), stdev_mg(ig, blitz::Range::all()), dims);
        }
    }

    // Fine flux tallies (TL)
    {
        ArrayB2 flux_mg(xs_mesh_.n_group(), mesh_.n_reg(MeshTreatment::TRUE));
        ArrayB2 stdev_mg(xs_mesh_.n_group(), mesh_.n_reg(MeshTreatment::TRUE));

        auto g = node.create_group("fsr_flux");
        for (int ig = 0; ig < (int)fine_flux_tally_.size(); ig++) {
            auto flux_result = fine_flux_tally_[ig].get();
            int ipin         = 0;
            for (const auto &v : flux_result) {
                flux_mg(ig, ipin)  = v.first;
                stdev_mg(ig, ipin) = v.second;
                ipin++;
            }
        }

        for (int ig = 0; ig < (int)fine_flux_tally_.size(); ig++) {
            std::stringstream path;
            path << std::setfill('0') << std::setw(3) << ig + 1;

            g.write(path.str(), flux_mg(ig, blitz::Range::all()));
            path << "_stdev";
            g.write(path.str(), stdev_mg(ig, blitz::Range::all()));
        }
    }

    // Fine flux tallies (collision)
    {
        VecF flux_mg(mesh_.n_reg(MeshTreatment::TRUE));
        VecF stdev_mg(mesh_.n_reg(MeshTreatment::TRUE));

        auto g = node.create_group("fsr_flux_col");
        for (int ig = 0; ig < (int)fine_flux_col_tally_.size(); ig++) {
            std::stringstream path;
            path << std::setfill('0') << std::setw(3) << ig + 1;
            auto flux_result = fine_flux_col_tally_[ig].get();
            int ireg         = 0;
            for (const auto &v : flux_result) {
                flux_mg[ireg]  = v.first;
                stdev_mg[ireg] = v.second;
                ireg++;
            }

            g.write(path.str(), flux_mg);
            path << "_stdev";
            g.write(path.str(), stdev_mg);
        }
    }

    // Pin powers
    {
        VecF pin_power;
        VecF stdev;
        pin_power.reserve(mesh_.n_pin());
        stdev.reserve(mesh_.n_pin());

        auto result = pin_power_tally_.get();
        for (const auto &v : result) {
            pin_power.push_back(v.first);
            stdev.push_back(v.second);
        }

        Normalize(pin_power.begin(), pin_power.end());

        node.write("pin_power", pin_power, dims);
        node.write("pin_power_stdev", stdev, dims);
    }

    return;
}
} // namespace mc
} // namespace mocc
