#pragma once

#include <iostream>
#include <vector>

#include "core_mesh.hpp"
#include "output_interface.hpp"

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
    // Using a structure of arrays would be more useful.
    class XSMeshRegion {
    friend class XSMesh;
    public:
        XSMeshRegion(){ }
        XSMeshRegion( const VecI& fsrs, const VecF& xstr, const VecF& xsnf,
                const VecF& xsch, const VecF& xsf, const ScatteringMatrix& xssc ):
            reg_(fsrs),
            xsmactr_(xstr),
            xsmacnf_(xsnf),
            xsmackf_(xsf),
            xsmacch_(xsch),
            xsmacsc_(xssc)
        {
            return;
        }

        size_t n_group() const {
            return xsmacsc_.n_group();
        }

        const real_t* xsmactr() const {
            return xsmactr_.data();
        }

        const real_t* xsmacnf() const {
            return xsmacnf_.data();
        }

        const real_t* xsmackf() const {
            return xsmackf_.data();
        }

        const real_t* xsmacch() const {
            return xsmacch_.data();
        }

        const ScatteringMatrix& xsmacsc() const {
            return xsmacsc_;
        }

        /**
         * Return a vector containing ALL of the FSRs that are filled with this
         * material.
         */
        const VecI& reg() const {
            return reg_;
        }

        friend std::ostream& operator<<(std::ostream& os,
                const XSMeshRegion &xsr );

    private:
        // List of FSR indices that use this XS mesh region
        VecI reg_;

        // Actual group constants for this XS mesh region
        VecF xsmactr_;
        VecF xsmacnf_;
        VecF xsmackf_;
        VecF xsmacch_;

        // Scattering matrix
        ScatteringMatrix xsmacsc_;
    };

    class XSMesh: public HasOutput {
    public:
        // Default constructor does almost nothing, and lets some other code
        // tell it what to do
        XSMesh() {}

        // XSMesh provides its own facility to initialize itself from a \ref
        // CoreMesh
        XSMesh( const CoreMesh& mesh );

        // Return the number of energy groups
        size_t n_group() const {
            return ng_;
        }

        // Iterators to the underlying vector
        const std::vector<XSMeshRegion>::const_iterator begin() const {
            return regions_.cbegin();
        }

        const std::vector<XSMeshRegion>::const_iterator end() const {
            return regions_.cend();
        }

        const XSMeshRegion& operator[]( size_t i ) const {
            return regions_[i];
        }

        size_t size() const {
            return regions_.size();
        }

        const VecF& eubounds() const {
            return eubounds_;
        }

        virtual void output( H5::CommonFG *file ) const {
            // Not really implementing for the general XS Mesh type.
            assert(false);
        }

    protected:
        size_t ng_;

        // Vector of xs mesh regions
        std::vector<XSMeshRegion> regions_;

        // Energy group upper bounds
        VecF eubounds_;
    };

    typedef std::shared_ptr<XSMesh> SP_XSMesh_t;
}
