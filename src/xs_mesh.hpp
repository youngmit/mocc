#pragma once

#include <vector>

#include "core_mesh.hpp"

// TODO:
//  - make a struct/class for an XS Mesh "region"
//  - make begin() and end() iterators
//  - the region class/struct can store data however it wants, but the interface
//  should be such that i can change the underlying structure later. Options for
//  that:
//      + store XS data in independent arrays for each region. This is easy to
//      implement and hard to screw up, but wont be the most efficient, since
//      the data will be all over the heap and carry around a bunch of overhead
//      + store all of the XS data on the xsmesh itself, and alias the proper
//      portion of the data in that blob. data will be super local, less
//      overhead, but a pain to implement.


namespace mocc {

    // For now im using the lazy implementation of the xsmesh region class.
    class XSMeshRegion {
    friend class XSMesh;
    public:
        XSMeshRegion( unsigned int ng, const VecI& fsrs ) {
            xsmactr_ = VecF(ng);
            xsmacnf_ = VecF(ng);
            xsmackf_ = VecF(ng);
            xsmacch_ = VecF(ng);
            reg_ = fsrs;
        }

        const float_t* xsmactr() const {
            return xsmactr_.data();
        }

        const float_t* xsmacnf() const {
            return xsmacnf_.data();
        }

        const float_t* xsmackf() const {
            return xsmackf_.data();
        }

        const float_t* xsmacch() const {
            return xsmacch_.data();
        }

        const VecI& reg() const {
            return reg_;
        }

    private:
        // List of FSR indices that use this XS mesh region
        VecI reg_;

        // Actual group constants for this XS mesh region
        VecF xsmactr_;
        VecF xsmacnf_;
        VecF xsmackf_;
        VecF xsmacch_;
    };

    class XSMesh {
    public:
        // Default constructor does almost nothing, and lets some other code
        // tell it what to do
        XSMesh() {}

        // XSMesh provides its own facility to initialize itself from a CoreMesh
        XSMesh( const CoreMesh& mesh );

        // Return the number of energy groups
        unsigned int n_grp() const {
            return ng_;
        }
        
        // Iterators to the underlying vector
        const std::vector<XSMeshRegion>::const_iterator begin() const {
            return regions_.cbegin();
        }

        const std::vector<XSMeshRegion>::const_iterator end() const {
            return regions_.cend(); 
        }

    private:
        // list of mesh regions corresponding to each XSMesh region
        unsigned int ng_;

        // Vector of xs mesh regions
        std::vector<XSMeshRegion> regions_;
    };
}
