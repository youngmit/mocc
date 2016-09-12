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

#include <iosfwd>
#include <vector>
#include "util/global_config.hpp"
#include "geometry/direction.hpp"
#include "geometry/geom.hpp"
#include "lattice.hpp"

namespace mocc {
class Plane {
public:
    Plane(const std::vector<const Lattice *> &lattices, size_t nx, size_t ny);

    const Lattice &at(size_t ix, size_t iy) const
    {
        assert(ix >= 0);
        assert(iy >= 0);
        assert(ix < nx_);
        assert(iy < ny_);
        return *(lattices_[ix + nx_ * iy]);
    }

    auto begin()
    {
        return lattices_.begin();
    }

    auto end()
    {
        return lattices_.end();
    }

    auto begin() const
    {
        return lattices_.cbegin();
    }

    auto end() const
    {
        return lattices_.cend();
    }

    /**
     * \brief Given a \ref Point2 in core-local coordinates, return a const
     * pointer to the corresponding \ref PinMesh.
     *
     * \param[in,out] p a Point2 in core-local coordinates. Will be modified
     * (see below).
     * \param[in,out] first_reg the first FSR index of the Plane. Will be
     * updated to the first region of the \ref Pin in which the pin resides.
     * \param[in] dir optional \ref Direction to use to disambiguate when
     * the \p p lies directly on a border.
     */
    const PinMesh *get_pinmesh(Point2 &p, int &first_reg,
                               Direction dir = Direction()) const;

    /**
     * \brief Return a const pointer to the \ref PinMesh that is at the
     * passed \ref Position.
     */
    const PinMesh *get_pinmesh(Position pos) const;

    /**
     * \brief Return the number of solution mesh regions in the \ref Plane
     */
    size_t n_reg() const
    {
        return n_reg_;
    }

    unsigned nx_pin() const
    {
        return nx_pin_;
    }

    unsigned ny_pin() const
    {
        return ny_pin_;
    }

    /**
     * \brief Return the total number of pins in all of the \ref Lattice objects
     * in this \ref Plane
     */
    size_t n_pin() const
    {
        return n_pin_;
    }

    /**
     * Return the number of XS Mesh regions in the Plane
     */
    size_t n_xsreg() const
    {
        return n_xsreg_;
    }

    /**
     * \brief Return a vector containing the FSR areas
     */
    VecF areas() const
    {
        VecF areas;
        for (auto &lat : lattices_) {
            for (auto &pin : *lat) {
                areas.insert(areas.end(), pin->areas().begin(),
                             pin->areas().end());
            }
        }

        return areas;
    }

    /**
     * \brief Return the position of a pin, given its index.
     *
     * The passed index will be cast into the valid range for the number of pins
     * in the plane by modulo. This allows pin indices larger than the number of
     * pins to still work, which is convenient in multi-plane situations where
     * the plane dimensions are identical.
     */
    Position pin_position(size_t ipin) const;

    /**
     * \brief Return the number of \ref Pin s in this \ref Plane marked as
     * fuel.
     */
    int n_fuel() const
    {
        return n_fuel_;
    }

    /**
     * \brief Return whether or not another \ref Plane is geometrically
     * identical to this \ref Plane
     *
     * This only checks PinMesh IDs; if there are multiple pin meshes with
     * different IDs but the same actual mesh structure, this method will
     * consider them to be different.
     */
    bool geometrically_equivalent(const Plane &other) const;

private:
    /**
     * Local list of \ref Lattice pointers
     */
    std::vector<const Lattice *> lattices_;

    /**
     * Number of lattices in the x direction
     */
    unsigned nx_;
    /**
     * Number of lattices in the y direction
     */
    unsigned ny_;

    unsigned nx_pin_;
    unsigned ny_pin_;

    size_t n_reg_;
    size_t n_xsreg_;

    /**
     * Locations of \ref Lattice interfaces
     */
    VecF hx_;
    VecF hy_;

    /**
     * List of the starting FSR index for each \ref Lattice in the plane
     */
    VecI first_reg_lattice_;

    int n_fuel_;

    int n_pin_;
};

/**
 * \brief Representation of a logical collection of Plane objects
 *
 * This struct contains useful data and references to data for the purpose of
 * interacting with the core mesh at the macroplane level. A collection of these
 * is intended to be very useful for the purpose of iterating over homogenized
 * planes in the mesh, and their internal pins.
 *
 * The \ref Plane reference contained within a \ref MacroPlane does not
 * necessarily contain the material data which should correspond to the \ref
 * MacroPlane, and should therefore not be used to obtain such data. Instead,
 * use the \ref begin() and \ref end() iterators to access the actual \ref Pin
 * objects and their \ref Material data.
 *
 * \todo The above is pretty dangerous. Might be good to provide a PlaneMesh
 * class which only stores PinMesh references, rather than Pin references.
 */
struct MacroPlane {
public:
    typedef std::vector<const Pin *>::const_iterator PinIter_t;
    MacroPlane(Plane const *plane, int iz_min, int iz_max, real_t height,
               PinIter_t begin, PinIter_t end)
        : plane(plane),
          iz_min(iz_min),
          iz_max(iz_max),
          height(height),
          begin_(begin),
          end_(end)
    {
        return;
    }

    PinIter_t begin() const
    {
        return begin_;
    }
    PinIter_t end() const
    {
        return end_;
    }
    const Pin *back() const
    {
        return *(end_ - 1);
    }
    unsigned size() const
    {
        return end_ - begin_;
    }
    unsigned size_3d() const
    {
        return (end_ - begin_) * (iz_max - iz_min + 1);
    }

    friend std::ostream &operator<<(std::ostream &os, const MacroPlane &mp);

    Plane const *plane;
    int iz_min;
    int iz_max;
    real_t height;

private:
    const PinIter_t begin_;
    const PinIter_t end_;
};
}
