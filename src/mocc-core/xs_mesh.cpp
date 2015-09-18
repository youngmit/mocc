#include "xs_mesh.hpp"

#include <iostream>
#include <memory>

#include "global_config.hpp"

using std::cout;
using std::endl;

namespace mocc {
    XSMesh::XSMesh( const CoreMesh& mesh ) {
        const MaterialLib& mat_lib = mesh.mat_lib();

        // Assume the same number of groups as the source material library
        ng_ = mat_lib.n_group();

        // Get energy group bounds
        eubounds_ = mat_lib.g_bounds();

        // loop over all of the pins in the CoreMesh and set up the XSMesh
        // regions
        std::vector<VecI> fsrs(mat_lib.n_materials());
        int ireg = 0;
        for( auto &pini: mesh ) {
            const PinMesh& pm = pini->mesh();
            const VecI& mat_ids = pini->mat_ids();
            int ixsreg = 0;
            for( auto &mat_id: mat_ids ) {
                int mat_index = mat_lib.get_index_by_id(mat_id);
                // post-increment pushes value, then increments
                for( size_t reg=0; reg<pm.n_fsrs(ixsreg); reg++) {
                    fsrs[mat_index].push_back(ireg++);
                }
                ixsreg++;
            }
        }

        // Calculate the necessary cross sections and store them in the
        // XSMesh-local arrays
        int imat = 0;
        for( auto &mat: mat_lib ) {
            // This is a bit of a hack until i use Eigen vectors or something
            // homebrew for storing these cross sections. need to add total
            // scattering to absorption to get xstr
            regions_.emplace_back( fsrs[imat], 
                    mat.xstr(),
                    mat.xsnf(), 
                    mat.xsch(), 
                    mat.xskf(), 
                    mat.xssc() );
            imat++;
        }

        return;
    }
}
