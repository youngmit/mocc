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
            // mi is a std::pair<int, Material>&
            matv.push_back(mi.first);
        }

        // Make a temporary map from the material ID as specified in the input
        // file/material library to its index in the XS mesh
        std::map<int, int> mat_map;
        int ixsreg = 0;
        for ( auto &i: matv ){
            mat_map[i] = ixsreg++;
        }

        // loop over all of the pins in the CoreMesh and set up the XSMesh
        std::vector<VecI> fsrs(mat_lib.n_materials());
        int ireg = 0;
        for( auto pini=mesh.begin_pin(); pini!=mesh.end_pin(); ++pini ) {
            const PinMesh& pm = (*pini)->mesh();
            const VecI& mat_ids = (*pini)->mat_ids();
            int ixsreg = 0;
            for( auto &mat: mat_ids ) {
                // post-increment pushes value, then increments
                for( int reg=0; reg<pm.n_fsrs(ixsreg); reg++) {
                    fsrs[mat].push_back(ireg++);
                }
            }
        }

        // lets see how we did
        for( auto &row: fsrs ) {
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
            regions_.emplace_back( ng_, fsrs[imat] );
            for( unsigned int ig=0; ig<ng_; ig++ ) {
                regions_.back().xsmacnf_[ig] = mati.xsnf()[ig];
                regions_.back().xsmactr_[ig] = mati.xsab()[ig]+mati.xssc().out(ig);
                regions_.back().xsmacch_[ig] = mati.xsch()[ig];
                regions_.back().xsmackf_[ig] = mati.xsf()[ig];
            }
        }

        // Again, lets see how we did
        //
        //for (auto &xsr: *this) {
        //    cout << "xsmesh region" << endl;
        //    cout << "ig      \t xstr    \t xsnf    \t xskf    \t xsch" << endl;
        //    for( int ig=0; ig<ng_; ig++ ) {
        //        cout << ig << "\t" 
        //             << xsr.xsmactr()[ig] << "\t"
        //             << xsr.xsmacnf()[ig] << "\t"
        //             << xsr.xsmackf()[ig] << "\t"
        //             << xsr.xsmacch()[ig] << "\t"
        //             << endl;
        //    }
        //}

        return;
    }
}
