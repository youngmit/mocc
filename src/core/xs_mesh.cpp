#include "xs_mesh.hpp"

#include <iostream>
#include <memory>
#include <map>

#include "blitz_typedefs.hpp"
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
        std::map<int, VecI> fsrs;
        int ireg = 0;
        for( auto &pini: mesh ) {
            const PinMesh& pm = pini->mesh();
            const VecI& mat_ids = pini->mat_ids();
            int ixsreg = 0;
            for( auto &mat_id: mat_ids ) {
                // post-increment pushes value, then increments
                for( size_t reg=0; reg<pm.n_fsrs(ixsreg); reg++) {
                    fsrs[mat_id].push_back(ireg++);
                }
                ixsreg++;
            }
        }

        int n_xsreg = fsrs.size();

        // Calculate the necessary cross sections and store them in the
        // XSMesh-local arrays
        xstr_.resize(n_xsreg, ng_);
        xsnf_.resize(n_xsreg, ng_);
        xsch_.resize(n_xsreg, ng_);
        xskf_.resize(n_xsreg, ng_);
        xsrm_.resize(n_xsreg, ng_);
        assert(xstr_(0, blitz::Range::all()).isStorageContiguous());

        VecI mat_ids(n_xsreg);
        int imat = 0;
        for( const auto &mat_pair: fsrs ) {
            mat_ids[imat] = mat_pair.first;
            const auto &mat = mat_lib[mat_pair.first];
            xstr_(imat, blitz::Range::all()) = mat.xstr();
            xsnf_(imat, blitz::Range::all()) = mat.xsnf();
            xsch_(imat, blitz::Range::all()) = mat.xsch();
            xskf_(imat, blitz::Range::all()) = mat.xskf();

            for( int ig=0; ig<(int)ng_; ig++ ) {
                xsrm_(imat, ig) = xstr_(imat, ig) - mat.xssc().self_scat(ig);
            }
            imat++;
        }
        
        // Preallocate space for the regions. Saves on lots of copies for large
        // xsmeshes.
        regions_.reserve(fsrs.size());
        imat = 0;
        for( auto &reg_pair: fsrs ) {
            const auto &mat = mat_lib[mat_ids[imat]];
            regions_.emplace_back( reg_pair.second,
                    &xstr_(imat, 0),
                    &xsnf_(imat, 0),
                    &xsch_(imat, 0),
                    &xskf_(imat, 0),
                    &xsrm_(imat, 0),
                    mat.xssc() );
            imat++;
        }

        return;
    }

    std::ostream& operator<<( std::ostream& os, const XSMeshRegion &xsr ) {
        int ng = xsr.xsmacsc_.n_group();
        os << "Transport: " << std::endl;
        for( int ig=0; ig<ng; ig++ ) {
            os << xsr.xsmactr_[ig];
        }
        os << std::endl;

        os << "nu-fission: " << std::endl;
        for( int ig=0; ig<ng; ig++ ) {
            os << xsr.xsmacnf_[ig];
        }
        os << std::endl;

        os << "chi: " << std::endl;
        for( int ig=0; ig<ng; ig++ ) {
            os << xsr.xsmacch_[ig];
        }
        os << std::endl;

        os << "removal: " << std::endl;
        for( int ig=0; ig<ng; ig++ ) {
            os << xsr.xsmacrm_[ig];
        }
        os << std::endl;

        os << "Scattering matrix:" << std::endl;
        os << xsr.xsmacsc_ << std::endl;

        return os;
    }
}
