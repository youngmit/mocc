#pragma once

#include <iostream>
#include <vector>

#include "blitz_typedefs.hpp"
#include "core_mesh.hpp"
#include "global_config.hpp"
#include "output_interface.hpp"

namespace mocc {
    // For now im using the lazy implementation of the xsmesh region class.
    // Using a structure of arrays would be more useful.
    class XSMeshRegion {
    friend class XSMesh;
    public:
        XSMeshRegion(){ }
        XSMeshRegion( const VecI &fsrs, real_t *xstr, 
                                        real_t *xsnf,
                                        real_t *xsch,
                                        real_t *xsf,
                                        real_t *xsrm,
                                        const ScatteringMatrix& xssc ):
            reg_(fsrs),
            xsmactr_(xstr),
            xsmacnf_(xsnf),
            xsmackf_(xsf),
            xsmacch_(xsch),
            xsmacrm_(xsrm),
            xsmacsc_(xssc)
        {
            return;
        }

        size_t n_group() const {
            return xsmacsc_.n_group();
        }

        const real_t* xsmactr() const {
            return xsmactr_;
        }

        const real_t* xsmacnf() const {
            return xsmacnf_;
        }

        const real_t* xsmackf() const {
            return xsmackf_;
        }

        const real_t* xsmacch() const {
            return xsmacch_;
        }

        const real_t* xsmacrm() const {
            return xsmacrm_;
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

        /**
         * \brief Update the cross sections referenced by the \ref XSMeshRegion
         */
        void update( const VecF &xstr, const VecF &xsnf, const VecF &xsch, 
                const VecF &xskf, const ScatteringMatrix &xssc )
        {
            for( int ig=0; ig<(int)xssc.n_group(); ig++ ) {
                xsmactr_[ig] = xstr[ig];
                xsmacnf_[ig] = xsnf[ig];
                xsmacch_[ig] = xsch[ig];
                xsmackf_[ig] = xskf[ig];
                xsmacrm_[ig] = xstr[ig] - xssc.self_scat(ig);
            }
            xsmacsc_ = xssc;
            return;
        }

        friend std::ostream& operator<<(std::ostream& os,
                const XSMeshRegion &xsr );

    private:
        // List of FSR indices that use this XS mesh region
        VecI reg_;

        // Actual group constants for this XS mesh region
        real_t *xsmactr_;
        real_t *xsmacnf_;
        real_t *xsmackf_;
        real_t *xsmacch_;
        real_t *xsmacrm_;

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

        // Actual cross-section data
        ArrayB2 xstr_;
        ArrayB2 xsnf_;
        ArrayB2 xsch_;
        ArrayB2 xskf_;
        ArrayB2 xsrm_;

        // Energy group upper bounds
        VecF eubounds_;
    };

    typedef std::shared_ptr<XSMesh> SP_XSMesh_t;
}
