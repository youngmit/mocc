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

#include <memory>

#include "util/blitz_typedefs.hpp"
#include "util/h5file.hpp"
#include "util/pugifwd.hpp"
#include "core/eigen_interface.hpp"
#include "core/xs_mesh.hpp"

namespace mocc {
/**
 * This class provides a homogenized version of an \ref XSMesh. At the time
 * of its conception, it was mostly intended for use with the Sn sweeper and
 * for future application to CMFD acceleration, so a lot of the code and
 * documentation make that assumption. If this is generalized much in the
 * future, make sure to update the docs.
 */
class XSMeshHomogenized : public XSMesh {
public:
    /**
     * Construct a homogenized \ref XSMesh from only a \ref CoreMesh object.
     * Following construction, the cross sections in the mesh will be
     * volume-weighted.
     */
    XSMeshHomogenized(const CoreMesh &mesh);

    /**
     * Construct a homogenized \ref XSMesh from a \ref CoreMesh and external
     * data. Following construction, the cross sections will be those
     * specified by the auxiliary data files specified by the \p input XML.
     *
     * \param mesh a \ref CoreMesh with which the \ref XSMeshHomogenized
     * should conform.
     * \param input an XML node containing at least one \<data\> child
     */
    XSMeshHomogenized(const CoreMesh &mesh, const pugi::xml_node &input);

    /**
     * \breif Construct a homogenized \ref XSMesh as a map into another
     * homogenized \ref XSMesh
     *
     * This allows two consumers of homogenized cross sections to use them on
     * different meshes. The new XSMesh will refer to the xs data of the
     * original XSMesh, but its XSMeshRetgions may have different metadata.
     */
    XSMeshHomogenized(XSMeshHomogenized &other,
                      const std::vector<int> &regions);

    /**
     * \brief \copybrief XSMesh::update()
     *
     * Update homogenized cross sections using the internally-stored
     * reference to the scalar flux.
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
    void set_flux(const ArrayB2 &flux)
    {
        // Make the assumption that the provided flux is using a PLANE treatment
        assert(flux.extent(0) == (int)mesh_.n_reg(MeshTreatment::PLANE));
        flux_ = &flux;
    }

    /**
     * Generate output of important cross sections on the homogenized mesh
     */
    void output(H5Node &file) const;

private:
    const CoreMesh &mesh_;

    // Possibly-associated flux for homogenization.
    const ArrayB2 *flux_;

    /**
    * \brief Populate the passed XSMeshRegion with homogenized cross sections
    * from a pin cell. No flux wieghting is performed, only volume weighting.
    *
    * \param i the index of the \ref XSMeshRegion in this \ref XSMeshHomogenized
    * \param pin the \ref Pin to use for the homogenization
    * \param xsr a reference to the \ref XSMeshRegion in which to store the
    * resulting
    * homogenized cross sections
    */
    void homogenize_region(int i, const Pin &pin, XSMeshRegion &xsr) const;

    /**
    * \brief Return an XSMeshRegion containing homogenized cross sections
    * from a pin cell. Use the passed scalar flux to perform flux-volume
    * weighting.
    *
    * \param i the index of the xs mesh region in this xsmesh to compute cross
    * sections for
    * \param first_reg the region offset into the \c flux array to be used
    * for this particular pin.
    * \param pin the \ref Pin to homogenize cross sections for.
    * \param xsr A reference to the \ref XSMeshRegion into which the cross
    * sections should be homogenized.
    *
    * This routine performs a flux-weighted homogenization of the cross
    * sections in the passed \ref Pin object and returns an \ref
    * XSMeshRegion object containing the homogenized cross sections.
    */
    void homogenize_region_flux(int i, int first_reg, const Pin &pin,
                                XSMeshRegion &xsr) const;

    void read_data_single(const pugi::xml_node &data);

    void read_data_multi(const pugi::xml_node &input);
};

typedef std::shared_ptr<XSMeshHomogenized> SP_XSMeshHomogenized_t;
}
