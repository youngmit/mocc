#pragma once

#include <memory>

#include "blitz_typedefs.hpp"
#include "eigen_interface.hpp"
#include "h5file.hpp"
#include "xs_mesh.hpp"

namespace mocc {
    /**
     * This class provides a homogenized version of an \ref XSMesh. At the time
     * of its conception, it was mostly intended for use with the Sn sweeper and
     * for future application to CMFD acceleration, so a lot of the code and
     * documentation make that assumption. If this is generalized much in the
     * future, make sure to update the docs.
     */
    class XSMeshHomogenized: public XSMesh {
    public:
        XSMeshHomogenized( const CoreMesh& mesh );

        /**
         * Update homogenized cross sections using the passed flux array
         */
        void update();

        /**
         * \brief Associate the homogenized cross-section mesh with a flux array
         * for flux-volume weighting.
         *
         * This associates the internal \c flux_ pointer with a fine mesh flux
         * array for use in the \ref XSMeshHomogenized::update() method. If no
         * flux has been associated when \ref update() is called, the update is
         * skipped, preserving the volume-weighted cross sections, which are
         * calculated at construction time.
         */
        void set_flux( const ArrayB2 &flux ) {
            assert( flux.extent(0) == (int)mesh_.n_reg() );
            flux_ = &flux;
        }

        /**
         * Generate output of important cross sections on the homogenized mesh
         */
        void output( H5::CommonFG *file ) const;
    private:
        const CoreMesh& mesh_;

        // Possibly-associated flux for homogenization.
        const ArrayB2 *flux_;


        /**
        * \brief Return an XSMeshRegion containing homogenized cross sections
        * from a pin cell. No flux wieghting is performed, only volume
        * weighting.
        */
        XSMeshRegion homogenize_region( int i, const Pin& pin ) const;

        /**
        * \brief Return an XSMeshRegion containing homogenized cross sections
        * from a pin cell. Use the passed scalar flux to perform flux-volume
        * weighting.
        *
        * \param i the region in the Sn or coarse mesh that the region should
        * belong to. Since it is assumed that there is a one-to-one mapping from
        * the xs mesh to the Sn mesh, the vector of FSRs in the \ref
        * XSMeshRegion will only contain one element, populated with the value
        * of \c i.
        * \param first_reg the region offset into the \c flux array to be used
        * for this particular pin.
        * \param pin the pin to homogenize cross sections for.
        * \param flux the array containing the scalar flux to be used in the
        * homogenization
        *
        * This routine performs a flux-weighted homogenization of the cross
        * sections in the passed \ref Pin object and returns an \ref
        * XSMeshRegion object containing the homogenized cross sections.
        */
        XSMeshRegion homogenize_region_flux( int i, int first_reg,
                const Pin& pin ) const;
    };

    typedef std::shared_ptr<XSMeshHomogenized> SP_XSMeshHomogenized_t;
}
