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

#include <algorithm>
#include <map>
#include <memory>
#include "util/pugifwd.hpp"
#include "core/assembly.hpp"
#include "core/core.hpp"
#include "core/lattice.hpp"
#include "core/material_lib.hpp"
#include "core/mesh.hpp"
#include "core/pin.hpp"
#include "core/pin_mesh.hpp"
#include "core/plane.hpp"

namespace mocc {
/**
* The core mesh stores everything needed to represent the physical state
* of the system. \ref PinMesh's, \ref MaterialLib's, \ref Pin's, \ref
* Lattice's etc. The \ref CoreMesh is then used to perform complex
* operations like ray tracing, generation of coarse mesh, etc. A lot of the
* heavy lifting for input processing happens in the constructor, and the
* CoreMesh assumes ownership of a lot of the structures used to represent
* the system.
*
* Once the sturctures that are used to define the system in the input file
* are parsed in, the \ref CoreMesh determines the set of
* geometrically-unique planes, which reduces the memory cost to perform ray
* tracing considerably.
*/
class CoreMesh : public Mesh {
public:
    /** Construct a \ref CoreMesh from XML input. This routine is
     * responsible for parsing many of the tags in the XML document:
     * \<mesh\>, \<pin\>, \<material_lib\>, \<lattice\>, \<core\>
    */
    CoreMesh(const pugi::xml_node &input);

    ~CoreMesh();

    /**
     * \brief Return the number of regions to consider for the given
     * MeshTreatment
     */
    size_t n_reg(MeshTreatment treatment) const override final
    {
        int n_reg = 0;
        switch (treatment) {
        case MeshTreatment::TRUE:
            n_reg = n_reg_;
            break;
        case MeshTreatment::PLANE: {
            int iz = 0;
            for (const auto np : subplane_) {
                n_reg += planes_[unique_plane_ids_[iz]].n_reg();
                iz += np;
            }
        } break;
        case MeshTreatment::PIN:
            n_reg = this->n_pin();
            break;
        case MeshTreatment::PIN_PLANE:
            n_reg = nx_ * ny_ * subplane_.size();
        }
        return n_reg;
    }

    /**
    * Return the number of geometrically-unique planes
    */
    size_t n_unique_planes() const
    {
        return planes_.size();
    }

    /**
     * Return the number of groups in the material library.
     */
    size_t n_group() const
    {
        return mat_lib_.n_group();
    }

    /**
    * \brief Obtain a tuple containing the pin position and a reference to
    * the PinMesh that occupies the space at a point, within a given plane.
    *
    * \param[inout] p a \ref Point2 residing in the desired \ref Pin. The
    * location of the point will be updated to the location of the \ref
    * PinMesh origin. See the note below.
    * \param[in] iplane the unique plane index of the \ref Plane in which to
    * search for the \ref PinMesh.
    * \param[inout] first_reg the index offset to start with. This will be
    * incremented by thie index of the first mesh region in the resultant
    * \ref PinMesh. This parameter will typically be passed in as zero, or
    * as the index of the first region of the <tt>iz</tt>th plane. The
    * former would be useful for ray tracing purposes, where each
    * geometrically-unique plane is traced independently, and the FSR
    * indices are incremented at sweep time to the appropriate concrete
    * plane. The latter is useful when the actual index is desired.
    *
    * This routine provides a means by which to locate the \ref PinMesh
    * object that fills the space in which the passed Point resides. This is
    * useful during ray tracing to determine the geometry that needs to be
    * traced through each pin subdomain. See \ref moc::Ray::Ray() for and
    * example of its use.
    *
    * \note The \ref Point2 that is passed in will be modified! The new
    * location will be the \ref PinMesh origin, in the core-local (global)
    * coordinate
    * system. This is done because during the ray trace, the original vector
    * of points coming from \ref Mesh::trace() are in core-local
    * coordinates, while the \ref PinMesh::trace() routine needs its \ref
    * Point2 (s) to be defined in pin-local coordinates, since the \ref
    * PinMesh has no idea where it is in global space. By moving the point
    * to the origin of the \ref PinMesh (but in core-local coordinates), we
    * can simply offset the ray points by the new location of \p p to get
    * into pin-local coordinates.
    */
    const PinMeshTuple get_pinmesh(Point2 &p, size_t iplane,
                                   int &first_reg) const;

    struct LocationInfo {
    public:
        const PinMesh *pm;
        Point2 local_point;
        int reg_offset;
        Position pos;
        std::array<Point3, 2> pin_boundary;
    };

    LocationInfo get_location_info(Point3 p, Direction dir) const;

    /**
    * \brief Return a const reference to the \ref Plane located at the indicated
    * axial position.
    *
    * \param ip The unique index for the \ref Plane to fetch
    */
    const Plane &unique_plane(int ip) const
    {
        assert((0 <= ip) & (ip < (int)planes_.size()));
        return planes_[ip];
    }

    /**
     * \brief Return a const reference to the \ref Plane that fills the passed
     * axial region
     */
    const Plane &get_plane_by_axial_index(int iz) const
    {
        assert((0 <= iz) & (iz < nz_));
        return planes_[unique_plane_ids_[iz]];
    }

    /**
     * \brief Return a reference to the vector of \ref MacroPlane objects
     *
     * This is useful for iterating over the macroplane/subplane mesh in
     * sweepers that support or are aware of such structure.
     */
    const std::vector<MacroPlane> &macroplanes() const
    {
        return macroplanes_;
    }

    /**
     * \brief Return a reference to the vector of macroplane heights
     */
    const VecF &macroplane_heights() const
    {
        return macroplane_heights_;
    }

    /**
     * \brief Return whether the mesh is Pin modular.
     *
     * In this context, pin modular essentially means that all of the pin
     * meshes in the mesh have the same pitch.
     */
    bool is_pin_modular() const
    {
        bool x_good = std::all_of(dx_vec_.begin(), dx_vec_.end(), [&](auto &v) {
            return fp_equiv_ulp(v, dx_vec_[0]);
        });
        bool y_good = std::all_of(dy_vec_.begin(), dy_vec_.end(), [&](auto &v) {
            return fp_equiv_ulp(v, dy_vec_[0]);
        });
        return x_good && y_good;
    }

    /**
    * \brief Return a const iterator to the first \ref Pin in the \ref
    * CoreMesh.
    */
    std::vector<const Pin *>::const_iterator begin() const
    {
        return core_pins_.cbegin();
    }

    /**
    * \brief Return a const iterator past the last \ref Pin in the \ref
    * CoreMesh.
    */
    std::vector<const Pin *>::const_iterator end() const
    {
        return core_pins_.cend();
    }

    /**
    * \brief Return a const iterator to the first \ref Pin in the \ref
    * CoreMesh in plane iz.
    */
    std::vector<const Pin *>::const_iterator begin(int iz) const
    {
        return core_pins_.cbegin() + nx_ * ny_ * iz;
    }

    /**
    * \brief Return a const iterator past the last \ref Pin in the \ref
    * CoreMesh in plane iz.
    */
    std::vector<const Pin *>::const_iterator end(int iz) const
    {
        return core_pins_.cbegin() + nx_ * ny_ * (iz + 1);
    }

    /**
    * \brief Return a const reference to the material library.
    */
    const MaterialLib &mat_lib() const
    {
        return mat_lib_;
    }

    /**
    * \brief Return the index of the first FSR within the given plane.
    */
    size_t first_reg_plane(int iz) const
    {
        assert((0 <= iz) & (iz < nz_));
        return first_reg_plane_[iz];
    }

    /**
    * \brief Return a \ref Position, indicating the global position of a pin
    * in the core geometry.
    *
    * At some point it would be nifty to create a custom iterator class that
    * can return this, obviating the need to keep track of pin index when
    * iterating over the core.
    */
    Position pin_position(size_t ipin) const;

    /**
     * \brief Return the Core-local coordinates to the pin origin
     */
    Point2 pin_origin(size_t ipin) const;

    /**
     * \brief Return a const reference to the vector of plane IDs
     *
     * The returned array is of size nz_, each entry being the index of the
     * unique plane that occupies that axial region.
     */
    const VecI &unique_plane_ids() const
    {
        return unique_plane_ids_;
    }

    /**
     * \brief Return the unique plane index correponding to an actual mesh plane
     */
    int unique_plane_id(int iz) const
    {
        return unique_plane_ids_[iz];
    }

    /**
     * \brief Returns true if this \ref CoreMesh represents a 2-D problem.
     *
     * A 2-D problem in this context means that there is only one plane, and
     * that both axial boundary conditions are reflective.
     */
    bool is_2d() const
    {
        return (nz_ == 1) && (bc_[(int)Surface::TOP] == Boundary::REFLECT) &&
               (bc_[(int)Surface::BOTTOM] == Boundary::REFLECT);
    }

    /**
     * \brief Print the contents of the mesh to the passed stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const CoreMesh &mesh);

    /**
     * \brief Return the region index at the point specified
     *
     * \param p a Point3 containing the coordinates to look up
     *
     * This method tracks down through the mesh hierarchy to find the region
     * index for the specified point. It is not fast.
     */
    int region_at_point(Point3 p) const
    {
        int iz            = this->plane_index(p.z);
        const auto &plane = planes_[unique_plane_ids_[iz]];

        // Start with the first region for the plane
        int ireg = first_reg_plane_[iz];

        // Get a pointer to the appropriate pin mesh, set ireg to the
        // beginning of that instance of the mesh.
        Point2 p2d    = p.to_2d();
        const auto pm = plane.get_pinmesh(p2d, ireg);

        ireg += pm->find_reg(p2d);

        return ireg;
    }

    const VecF volumes(MeshTreatment treatment) const
    {
        VecF volumes;
        switch (treatment) {
        case MeshTreatment::TRUE: {
            volumes.reserve(this->n_reg(treatment));
            int ipin = 0;
            for (const auto &pin : *this) {
                real_t hz     = dz_vec_[ipin / (nx_ * ny_)];
                auto &pinmesh = pin->mesh();
                for (const auto &a : pinmesh.areas()) {
                    volumes.push_back(a * hz);
                }
                ipin++;
            }
            return volumes;
        }
        case MeshTreatment::PLANE: {
            VecF volumes;
            volumes.reserve(this->n_reg(treatment));

            int iz     = 0;
            int imacro = 0;
            for (const auto np : subplane_) {
                real_t hz = macroplane_heights_[imacro];

                auto stt = this->begin() + nx_ * ny_ * iz;
                auto stp = stt + nx_ * ny_;
                for (auto it = stt; it != stp; ++it) {
                    auto &pinmesh = (*it)->mesh();
                    for (const auto &a : pinmesh.areas()) {
                        volumes.push_back(a * hz);
                    }
                }

                iz += np;
                imacro++;
            }
            return volumes;
        }
        case MeshTreatment::PIN_PLANE:
            volumes.reserve(this->n_reg(treatment));
            for (const auto &mplane : macroplanes_) {
                for (int iy = 0; iy < ny_; iy++) {
                    for (int ix = 0; ix < nx_; ix++) {
                        volumes.push_back(dx_vec_[ix] * dy_vec_[iy] *
                                          mplane.height);
                    }
                }
            }

            return volumes;
        case MeshTreatment::PIN:
            // Fall through and return coarse volume. Putting the last return in
            // here causes warnings on some compilers
            break;
        }

        return this->coarse_volume();
    }

    const auto &pins() const
    {
        return pins_;
    }

    /**
     * \brief Return a reference to the subplane parameters
     */
    const VecI subplane() const
    {
        return subplane_;
    }

    const Core &core() const
    {
        return core_;
    }

    int n_fuel_2d() const
    {
        return n_fuel_2d_;
    }

private:
    // Map for storing pin mesh objects indexed by user-specified IDs
    std::map<int, UP_PinMesh_t> pin_meshes_;

    // The material library
    MaterialLib mat_lib_;

    // Map of actual pin objects, indexed by user-specified IDs
    std::map<int, UP_Pin_t> pins_;

    // Map of lattice objects
    std::map<int, UP_Lattice_t> lattices_;

    // Map of assembly objects
    std::map<int, UP_Assembly_t> assemblies_;

    // Vector of Plane instances. There should be one for each unique planar
    // geometry
    std::vector<Plane> planes_;

    // Vector of references to all pins in the core. This facilitates
    // iteration through the entire problem geometry in a linear fashion.
    // The pins are ordered in the same order that flat source regions are
    // indexed.
    std::vector<const Pin *> core_pins_;

    // Core object (essentially a 2D array of Assemblies)
    Core core_;

    // Subplane factors
    VecI subplane_;

    // The height of each macroplane
    VecF macroplane_heights_;

    std::vector<MacroPlane> macroplanes_;

    // Number of assemblies
    size_t nasy_;

    // Number of pins to consider "fuel." Useful for normalization
    int n_fuel_2d_;

    // List of geometrically-unique planes. Each entry in the list
    // corresponds to the unique plane index that is geometrically valid for
    // the actual plane.
    VecI unique_plane_ids_;

    // Plane index of the first occurance of each geometrically-unique plane
    VecI first_unique_;

    // Index of the first flat source region on each plane
    VecI first_reg_plane_;
};

typedef std::shared_ptr<CoreMesh> SP_CoreMesh_t;
typedef std::unique_ptr<CoreMesh> UP_CoreMesh_t;
}
