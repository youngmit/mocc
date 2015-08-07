#include "xs_mesh_homogenized.hpp"

namespace mocc {
    

    XSMeshHomogenized::XSMeshHomogenized( const CoreMesh& mesh ):
        mesh_( mesh )
    {
        // Set up the non-xs part of the xs mesh
        eubounds_ = mesh_.mat_lib().g_bounds();
        ng_ = eubounds_.size();
        
        int ipin = 0;
        int first_reg = 0;
        for( auto &pin: mesh_ ) {
            regions_.push_back(this->homogenize_region( ipin, first_reg, *pin ));
            ipin++;
            first_reg += pin->n_reg();
        }

    }

    XSMeshRegion XSMeshHomogenized::homogenize_region( int i, int first_reg,
            const Pin& pin) const {
        VecI fsrs( 1, i );
        VecF xstr( ng_, 0.0 );
        VecF xsnf( ng_, 0.0 );
        VecF xskf( ng_, 0.0 );
        VecF xsch( ng_, 0.0 );

        std::vector<VecF> scat( ng_, VecF(ng_, 0.0) );

        auto mat_lib = mesh_.mat_lib();
        auto &pin_mesh = pin.mesh();
        auto vols = pin_mesh.vol();

        for( int ig=0; ig<ng_; ig++ ) {
            int ireg = 0;
            int ixsreg = 0;
            for( auto &mat_id: pin.mat_ids() ) {
                auto mat = mat_lib.get_material_by_id(mat_id);
                for( int i=0; i<pin_mesh.n_fsrs(ixsreg); i++ ) {
                    xstr[ig] += vols[ireg] * mat.xsab()[ig];
                    ireg++;
                }
                ixsreg++;
            }
        }

        ScatMat scat_mat(scat);

        return XSMeshRegion( fsrs, xstr, xsnf, xsch, xskf, scat_mat );
    }
}
