#include "xs_mesh.hpp"

#include <iostream>

#include "global_config.hpp"

using std::cout;
using std::endl;

namespace mocc {
    
    XSMesh::XSMesh( const CoreMesh& mesh ) {
        const MaterialLib& mat_lib = mesh.mat_lib();

        // Assume the same number of groups as the source material library
        ng_ = mat_lib.n_grp();

        // Make a list of the material ids, to get back into a dense index space
        VecI matv;
        for( auto &mi: mat_lib.materials()) {
            matv.push_back(mi.first);
        }

        // loop over all of the pins in the CoreMesh and set up the XSMesh
        reg_ = std::vector<VecI>(mat_lib.n_materials());
        int ireg = 0;
        for( auto pini=mesh.begin_pin(); pini!=mesh.end_pin(); ++pini ) {
            const PinMesh& pm = (*pini)->mesh();
            const VecI& mat_ids = (*pini)->mat_ids();
            int ixsreg = 0;
            for( auto &mat: mat_ids ) {
                // post-increment pushes value, then increments
                for( int reg=0; reg<pm.n_fsrs(ixsreg); reg++) {
                    reg_[mat].push_back(ireg++);
                }
            }
        }

        // lets see how we did
        for( auto &row: reg_ ) {
            cout << "xs mesh region fsrs: ";
            for( auto &i: row ) {
                cout << i << " ";
            }
            cout << endl;
        }

        // Calculate the necessary cross sections and store them in the
        // XSMesh-local arrays
        for( unsigned int imat=0; imat<mat_lib.n_materials(); imat++ ) {
            unsigned int lib_key = matv[imat];
            const Material& mati = *mat_lib.materials().at(lib_key);
            for( unsigned int ig=0; ig<ng_; ig++ ) {
                xsmacnf_.push_back(mati.xsnf()[ig]);
                xsmactr_.push_back(mati.xsab()[ig]+mati.xssc().out(ig));
                xsmacch_.push_back(mati.xsch()[ig]);
                xsmackf_.push_back(mati.xsf()[ig]);
            }
        }

        return;
    }
}
