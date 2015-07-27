#pragma once

#include <map>
#include <memory>

#include "pugixml.hpp"

#include "pin_mesh.hpp"
#include "pin.hpp"
#include "material_lib.hpp"
#include "lattice.hpp"
#include "assembly.hpp"
#include "core.hpp"
#include "plane.hpp"

namespace mocc {

    // The core mesh stores everything needed to represent the physical state
    // of the system. Pin meshes, material library, actual pin types, lattices
    // etc. The CoreMesh is then used to perform complex operations like ray
    // tracing, generation of coarse mesh, etc. A lot of the heavy lifting for
    // input processing happens in the constructor, and the CoreMesh assumes 
    // ownership of a lot of the structures used to represent the system.
    class CoreMesh {
    public:
        CoreMesh() {
            // do nothing
            return;
        }

        // Construct a CoreMesh from XML input. This routine is responsible for
        // parsing many of the tags in the XML document: <mesh>, <pin>,
        // <material_lib>, <lattice>, <core>
        CoreMesh( pugi::xml_node &input );

        ~CoreMesh();

        float_t hx() const {
            return hx_;
        }

        float_t hy() const {
            return hy_;
        }

        unsigned int nx() const {
            return nx_;
        }

        unsigned int ny() const {
            return ny_;
        }

        unsigned int nz() const {
            return nz_;
        }

        int n_unique_planes() const {
            return first_unique_.size();
        }

        float_t n_reg() const {
            return n_reg_;
        }
        
        // Given a vector containing two points (Which should be on the boundary
        // of the core mesh), insert points corresponding to intersections of
        // the line formed by those points. The points are added to the vector
        // itself.
        void trace( std::vector<Point2> &p ) const;           

        // return a reference to the PinMesh that occupies the space at a point,
        // within a given plane.
        const PinMesh* get_pinmesh( Point2 &p, unsigned int iz, 
                int &first_reg) const;

        const Plane& plane( unsigned int iz ) const {
            assert( (0 <= iz) & (iz < planes_.size()) );
            return planes_[iz];
        }

        std::vector<const Pin*>::const_iterator begin_pin() const {
            return core_pins_.cbegin();
        }

        std::vector<const Pin*>::const_iterator end_pin() const {
            return core_pins_.cend();
        }
        
        // Return a reference to the material library
        const MaterialLib& mat_lib() const {
            return mat_lib_;
        }

        // Return the first FSR index for a given plane
        unsigned int first_reg_plane( unsigned int iz ) const {
            assert( (0 <= iz ) & ( iz<nz_ ) );
            return first_reg_plane_[iz];
        }

        // Return the core boundary conditions
        std::vector<Boundary> boundary() const {
            return core_.boundary();
        }

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

        // Total core size in the x dimension
        float_t hx_;

        // Total core size in the y dimension
        float_t hy_;

        // Total core size in the z dimension
        float_t hz_;

        // List of pin boundaries in the x dimension
        VecF x_vec_;

        // List of pin boundaries in the y dimension
        VecF y_vec_;

        // Numbers of pins/planes in each dimension
        unsigned int nx_;
        unsigned int ny_;
        unsigned int nz_;

        // Number of assemblies
        unsigned int nasy_;

        // Total number of FSRs in the entire geometry
        unsigned int n_reg_;
        // Total number of XS regions in the entire geometry
        unsigned int n_xsreg_;

        // List of geometrically-unique planes. Each entry in the list
        // corresponds to the unique plane index that is geometrically valid for
        // the actual plane.
        VecI unique_plane_;

        // Plane index of the first occurance of each geometrically-unique plane
        VecI first_unique_;

        // Index of the first flat source region on each plane
        VecI first_reg_plane_;

        // Vector of Line objects, representing pin boundaries. This greatly
        // simplifies the ray trace.
        std::vector<Line> lines_;

        // Vector containing the flat source region index of the first region in
        // each plane.


    };

    typedef std::shared_ptr<CoreMesh> SP_CoreMesh_t;
    typedef std::unique_ptr<CoreMesh> UP_CoreMesh_t;
}
