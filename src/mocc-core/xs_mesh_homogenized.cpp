#include "xs_mesh_homogenized.hpp"

#include <iostream>

using std::cout;
using std::endl;

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
            // Use the lexicographically-ordered pin index as the xs mesh index.
            // This is to put the indexing in a way that works best for the Sn
            // sweeper as it is implemented now. This is really brittle, and
            // should be replaced with some sort of Sn Mesh object, which both
            // the XS Mesh and the Sn sweeper will use to handle indexing.
            int ireg = mesh.index_lex( mesh.pin_position(ipin) );
            regions_.push_back(this->homogenize_region( ireg, first_reg, 
                        *pin ));
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
        auto vols = pin_mesh.vols();

        for( unsigned int ig=0; ig<ng_; ig++ ) {
            int ireg = 0;
            int ixsreg = 0;
            for( auto &mat_id: pin.mat_ids() ) {
                auto mat = mat_lib.get_material_by_id(mat_id);
                const ScatRow& scat_row = mat.xssc().to(ig);
                int gmin = scat_row.min_g;
                int gmax = scat_row.max_g;
                for( int i=0; i<pin_mesh.n_fsrs(ixsreg); i++ ) {
                    xstr[ig] += vols[ireg] * mat.xstr()[ig];
                    xsnf[ig] += vols[ireg] * mat.xsnf()[ig];
                    xskf[ig] += vols[ireg] * mat.xskf()[ig];
                    xsch[ig] += vols[ireg] * mat.xsch()[ig];

                    for( int igg=gmin; igg<=gmax; igg++ ) {
                        scat[ig][igg] += scat_row.from[igg-gmin] * vols[ireg];
                    }
                    ireg++;
                }
                ixsreg++;
            }

            xstr[ig] /= pin.vol();
            xsnf[ig] /= pin.vol();
            xskf[ig] /= pin.vol();
            xsch[ig] /= pin.vol();

            for( auto &s: scat[ig] ) {
                s /= pin.vol();
            }

        }

        /// \todo i guess i never finished scattering matrix homogenization?
        //for( auto &row: scat ) {
        //    for( auto &c: row) {
        //    }
        //}

        ScatMat scat_mat(scat);

        return XSMeshRegion( fsrs, xstr, xsnf, xsch, xskf, scat_mat );
    }

    void update( const ArrayX &flux ){
        return;
    }
}
