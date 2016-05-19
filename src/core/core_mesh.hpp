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

#include "core/assembly.hpp"
#include "core/core.hpp"
#include "core/lattice.hpp"
#include "core/material_lib.hpp"
#include "core/mesh.hpp"
#include "core/plane.hpp"
#include "core/pin.hpp"
#include "core/pin_mesh.hpp"
#include "core/pugifwd.hpp"

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
    class CoreMesh: public Mesh {
    public:
        /** Construct a \ref CoreMesh from XML input. This routine is
         * responsible for parsing many of the tags in the XML document:
         * \<mesh\>, \<pin\>, \<material_lib\>, \<lattice\>, \<core\>
        */
        CoreMesh( const pugi::xml_node &input );

        ~CoreMesh();

        /**
        * Return the number of geometrically-unique planes
        */
        size_t n_unique_planes() const {
            return first_unique_.size();
        }

        /**
         * Return the number of groups in the material library.
         */
        size_t n_group() const {
            return mat_lib_.n_group();
        }

        /**
        * \brief Obtain a tuple containing the pin position and a reference to
        * the PinMesh that occupies the space at a point, within a given plane.
        *
        * \param[inout] p a \ref Point2 residing in the desired \ref Pin. The
        * location of the point will be updated to the location of the \ref
        * PinMesh origin. See the note below.
        * \param[in] iz the index of the \ref Plane in which to search for the
        * \ref PinMesh.
        * \param[inout] first_reg still need to accurately define
        * this param
        *
        * \todo document thie first_reg parameter. important!
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
        const PinMeshTuple get_pinmesh( Point2 &p, size_t iz,
                int &first_reg) const;

        /**
        * \brief Return a const reference to the indexed plane.
        *
        * These planes are not considered "unique;" it returns actual Plane that
        * fills the indezed axial region.
        */
        const Plane& plane( unsigned int iz ) const {
            assert( (0 <= iz) & (iz < planes_.size()) );
            return planes_[iz];
        }

        /**
         * \brief Return whether the mesh is Pin modular.
         *
         * In this context, pin modular essentially means that all of the pin
         * meshes in the mesh have the same pitch.
         */
        bool is_pin_modular() const {
            bool x_good = std::all_of( dx_vec_.begin(), dx_vec_.end(),
                    [&] (auto &v)
                    {
                        return fp_equiv_ulp(v, dx_vec_[0]);
                    }
                );
            bool y_good = std::all_of( dy_vec_.begin(), dy_vec_.end(),
                    [&] (auto &v)
                    {
                        return fp_equiv_ulp(v, dy_vec_[0]);
                    }
                );
            return x_good && y_good;
        }

        /**
        * \brief Return a const iterator to the first \ref Pin in the \ref
        * CoreMesh.
        */
        std::vector<const Pin*>::const_iterator begin_pin() const {
            return core_pins_.cbegin();
        }

        /**
        * \brief Return a const iterator past the last \ref Pin in the \ref
        * CoreMesh.
        */
        std::vector<const Pin*>::const_iterator end_pin() const {
            return core_pins_.cend();
        }

        /**
        * \brief Return a const iterator to the first \ref Pin in the \ref
        * CoreMesh.
        */
        std::vector<const Pin*>::const_iterator begin() const {
            return core_pins_.cbegin();
        }

        /**
        * \brief Return a const iterator past the last \ref Pin in the \ref
        * CoreMesh.
        */
        std::vector<const Pin*>::const_iterator end() const {
            return core_pins_.cend();
        }

        /**
        * \brief Return a const iterator to the first \ref Pin in the \ref
        * CoreMesh.
        */
        std::vector<const Pin*>::const_iterator begin( int iz ) const {
            return core_pins_.cbegin()+nx_*ny_*iz;
        }

        /**
        * \brief Return a const iterator past the last \ref Pin in the \ref
        * CoreMesh.
        */
        std::vector<const Pin*>::const_iterator end( int iz ) const {
            return core_pins_.cbegin()+nx_*ny_*(iz+1);
        }

        /**
        * \brief Return a const reference to the material library.
        */
        const MaterialLib& mat_lib() const {
            return mat_lib_;
        }

        /**
        * \brief Return the index of the first FSR within the given plane.
        */
        size_t first_reg_plane( int iz ) const {
            assert( (0 <= iz ) & ( iz < nz_ ) );
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
        Position pin_position( size_t ipin ) const;

        /**
         * \brief Return the Core-local coordinates to the pin origin
         */
        Point2 pin_origin( size_t ipin ) const;

        /**
         * \brief Return a const reference to the vector of plane IDs
         */
        const VecI& unique_planes() const {
            return unique_plane_;
        }

        /**
         * \brief Return the number of fuel pins in the plane.
         *
         * Practically speaking, this is the maximum number of pins in all of
         * the Planes of the core.
         */
        int n_fuel_2d() const {
            return n_fuel_2d_;
        }

        /**
         * \brief Returns true if this \ref CoreMesh represents a 2-D problem.
         *
         * A 2-D problem in this context means that there is only one plane, and
         * that both axial boundary conditions are reflective.
         */
        bool is_2d() const {
            return (nz_ == 1) && (bc_[Surface::TOP] == Boundary::REFLECT) &&
                (bc_[Surface::BOTTOM] == Boundary::REFLECT);
        }

        /**
         * \brief Print the contents of the mesh to the passed stream.
         */
        friend std::ostream& operator<<( std::ostream &os,
                const CoreMesh &mesh);

        /**
         * \brief Return the region index at the point specified
         *
         * \param p a Point3 containing the coordinates to look up
         *
         * This method tracks down through the mesh hierarchy to find the region
         * index for the specified point. It is not fast.
         */
        int region_at_point( Point3 p ) const {
            int plane_index = std::distance(z_vec_.begin(),
                    std::lower_bound(z_vec_.begin(), z_vec_.end(), p.z));
            const auto &plane = planes_[plane_index];

            // Start with the first region for the plane
            int ireg = first_reg_plane_[plane_index];

            // Get a pointer to the appropriate pin mesh, set ireg to the
            // beginning of that instance of the mesh.
            Point2 p2d = p.to_2d();
            const auto pm = plane.get_pinmesh( p2d, ireg );

            ireg += pm->find_reg(p2d);

            return ireg;
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
        std::vector<const Pin*> core_pins_;

        // Core object (essentially a 2D array of Assemblies)
        Core core_;

        // List of plane heights
        VecF hz_vec_;

        // Number of assemblies
        size_t nasy_;

        // List of geometrically-unique planes. Each entry in the list
        // corresponds to the unique plane index that is geometrically valid for
        // the actual plane.
        VecI unique_plane_;

        // Plane index of the first occurance of each geometrically-unique plane
        VecI first_unique_;

        // Index of the first flat source region on each plane
        VecI first_reg_plane_;

        int n_fuel_2d_;
    };

    typedef std::shared_ptr<CoreMesh> SP_CoreMesh_t;
    typedef std::unique_ptr<CoreMesh> UP_CoreMesh_t;
}
