#pragma once

#include <map>
#include <memory>

#include "pugixml.hpp"

#include "mesh.hpp"
#include "pin_mesh.hpp"
#include "pin.hpp"
#include "material_lib.hpp"
#include "lattice.hpp"
#include "assembly.hpp"
#include "core.hpp"
#include "plane.hpp"

namespace mocc {
    /**
    * The core mesh stores everything needed to represent the physical state
    * of the system. PinMesh's, MaterialLib's, Pin's, Lattice's
    * etc. The CoreMesh is then used to perform complex operations like ray
    * tracing, generation of coarse mesh, etc. A lot of the heavy lifting for
    * input processing happens in the constructor, and the CoreMesh assumes
    * ownership of a lot of the structures used to represent the system.
    *
    * Once the sturctures that are used to define the system in the input file
    * are parsed in, the CoreMesh determines the set of geometrically-unique
    * planes, which reduces the memory cost to perform ray tracing considerably.
    */
    class CoreMesh: public Mesh {
    public:
        /** Construct a CoreMesh from XML input. This routine is responsible for
         * parsing many of the tags in the XML document: \<mesh\>, \<pin\>,
         * \<material_lib\>, \<lattice\>, \<core\>
        */
        CoreMesh( pugi::xml_node &input );

        ~CoreMesh();

        /** 
        * Return the total core length along the x dimension
        */
        real_t hx() const {
            return hx_;
        }

        /** 
        * Return the total core length along the y dimension
        */
        real_t hy() const {
            return hy_;
        }

        /** 
        * Return the pin boundary locations along the x dimension
        */
        const VecF& pin_hx() const {
            return x_vec_;
        }

        /** 
        * Return the pin boundary locations along the y dimension
        */
        const VecF& pin_hy() const {
            return y_vec_;
        }

        /** 
        * Return a const vector of plane heights
        */
        const VecF& hz() const {
            return hz_vec_;
        }

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
        * \brief Obtain a tuple containing the pin position and a reference to the
        * PinMesh that occupies the space at a point, within a given plane.
        *
        * \param[in,out] p a Point residing in the desired Pin. The location of
        * the point will be updated to the location of the PinMesh origin. See
        * the note below.
        * \param iz[in] the index of the Plane in which to search for the PinMesh.
        * \param first_reg[in,out]
        *
        * This routine provides a means by which to locate the PinMesh object
        * that fills the space in which the passed Point resides. This is useful
        * during ray tracing to determine the geometry that needs to be traced
        * through each pin subdomain. See Ray::Ray() for and example of its use.
        *
        * \note The Point that is passed in will be modified! The new location
        * will be the PinMesh origin, in the core-local (global) coordinate
        * system. This is done because during the ray trace, the original vector
        * of points coming from CoreMesh::trace() are in core-local
        * coordinates, while the PinMesh::trace() routine needs its Point(s) to
        * be defined in pin-local coordinates, since the PinMesh has no idea
        * where it is in global space. By moving the point to the origin of the
        * PinMesh (but in core-local coordinates), we can simply offset the ray
        * points by the new location of \p p to get into pin-local coordinates.
        */
        const PinMeshTuple get_pinmesh( Point2 &p, unsigned int iz, 
                int &first_reg) const;

        /**
        * Return a const reference to the indexed plane. These planes are not
        * considered "unique;" it returns actual Plane that fills the indezed
        * axial region.
        */
        const Plane& plane( unsigned int iz ) const {
            assert( (0 <= iz) & (iz < planes_.size()) );
            return planes_[iz];
        }

        /**
        * Return a const iterator to the first Pin in the CoreMesh.
        */
        std::vector<const Pin*>::const_iterator begin_pin() const {
            return core_pins_.cbegin();
        }

        /**
        * Return a const iterator past the last Pin in the CoreMesh.
        */
        std::vector<const Pin*>::const_iterator end_pin() const {
            return core_pins_.cend();
        }

        /**
        * Return a const iterator to the first Pin in the CoreMesh.
        */
        std::vector<const Pin*>::const_iterator begin() const {
            return core_pins_.cbegin();
        }

        /**
        * Return a const iterator past the last Pin in the CoreMesh.
        */
        std::vector<const Pin*>::const_iterator end() const {
            return core_pins_.cend();
        }
        
        /**
        * Return a const reference to the material library.
        */
        const MaterialLib& mat_lib() const {
            return mat_lib_;
        }

        /**
        * Return the index of the first FSR within the given plane.
        */
        unsigned int first_reg_plane( unsigned int iz ) const {
            assert( (0 <= iz ) & ( iz<nz_ ) );
            return first_reg_plane_[iz];
        }

        /**
        * Return a vector containing the core boundary conditions.
        */
        std::vector<Boundary> boundary() const {
            return core_.boundary();
        }

        /**
        * Return the 1-D lexicographic index of a Position.
        *
        * It is often useful to flatten the index space of the entire geometry
        * for, usually to facilitate output to the HDF5 file. This accepts a
        * Position, which stores x, y, z indices for a Pin in the mesh and
        * returns a 1-D index. The 1-D index space is ordered lexicographically,
        * meaning in ascending x, then y, then z. Otherwise called "natural"
        * ordering.
        */
        unsigned int index_lex( Position pos ) const {
            return pos.x + pos.y*nx_ + pos.z*nx_*ny_;
        }
        
        /*
        * Return a Position, indicating the global position of a pin in the
        * core geometry.
        *
        * At some point it would be nifty to create a custom iterator class that
        * can return this, obviating the need to keep track of pin index when
        * iterating over the core.
        */
        Position pin_position( unsigned int ipin ) const;

    private:
        // Map for storing pin mesh objects indexed by user-specified IDs
        std::map<int, UP_PinMesh_t> pin_meshes_;

        // The material library
        MaterialLib mat_lib_;

        // Map of actual pin objects, indexed by user-specified IDs
        std::map<int, UP_Pin_t> pins_;

        // Map of lattice objects
        std::map<int, Lattice> lattices_;

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

        
    };

    typedef std::shared_ptr<CoreMesh> SP_CoreMesh_t;
    typedef std::unique_ptr<CoreMesh> UP_CoreMesh_t;
}
