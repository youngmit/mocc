#include "scattering_matrix.hpp"

#include <iomanip>
#include <iostream>

namespace mocc {
    ScatteringMatrix::ScatteringMatrix( const std::vector<VecF> &scat){

        // Imply ng_ from the size of the passed-in vectors
        ng_ = scat.size();
        out_ = VecF(ng_, 0.0);

        // Densify the scattering matrix
        for( int to=0; to<ng_; to++ ) {
            for( int from=0; from<ng_; from++ ) {
                if(scat[to][from] > 0.0){
                    scat_.push_back(scat[to][from]);
                }
                out_[from] += scat[to][from];
            }
        }

        // Determine group bounds for each row
        // TODO: This is some nasty jazz, might have originally written this
        // whilst under the influence of something... Come back and clean up
        // later.
        int pos = 0;
        int prevPos = 0;
        int min_g = 0;
        int max_g = 0;
        for( int to=0; to<ng_; to++ ) {
            bool found_min = false;
            bool found_max = false;
            for( int from=0; from<ng_; from++ ) {
                if( scat[to][from] > 0.0 ) {
                    if(!found_min){
                        found_min = true;
                        min_g = from;
                    }
                    pos++;
                }
                if( scat[to][from] == 0.0 && found_min ) {
                    found_max = true;
                    max_g = from-1;
                    break;
                }
            }
            if( !found_max ) {
                max_g = ng_-1;
            }
            rows_.push_back(ScatteringRow(min_g, max_g, &scat_[prevPos]));
            prevPos = pos;
        }
        return;
    }

    ScatteringMatrix::ScatteringMatrix(const ArrayB2 &scat){
        // Make sure that scat is square
        assert(scat.extent(0) == scat.extent(1));
        // Imply ng_ from the size of the passed-in vectors
        ng_ = scat.extent(0);
        out_ = VecF(ng_, 0.0);

        // Densify the scattering matrix
        for( int to=0; to<ng_; to++ ) {
            for( int from=0; from<ng_; from++ ) {
                if(scat(to, from) > 0.0){
                    scat_.push_back(scat(to, from));
                }
                out_[from] += scat(to, from);
            }
        }

        // Determine group bounds for each row
        // TODO: This is some nasty jazz, might have originally written this
        // whilst under the influence of something... Come back and clean up
        // later.
        int pos = 0;
        int prevPos = 0;
        int min_g = 0;
        int max_g = 0;
        for( int to=0; to<ng_; to++ ) {
            bool found_min = false;
            bool found_max = false;
            for( int from=0; from<ng_; from++ ) {
                if( scat(to, from) > 0.0 ) {
                    if(!found_min){
                        found_min = true;
                        min_g = from;
                    }
                    pos++;
                }
                if( scat(to, from) == 0.0 && found_min ) {
                    found_max = true;
                    max_g = from-1;
                    break;
                }
            }
            if( !found_max ) {
                max_g = ng_-1;
            }
            rows_.push_back(ScatteringRow(min_g, max_g, &scat_[prevPos]));
            prevPos = pos;
        }
        return;
    }

    std::ostream& operator<<( std::ostream& os,
            const ScatteringMatrix &scat_mat ) {
        for( auto &row: scat_mat.rows_ ) {
            int gmin = row.min_g;
            int gmax = row.max_g;
            for( int ig=0; ig<gmin; ig++ ) {
                os << std::setw(12) << 0.0;
            }

            for( int ig=gmin; ig<=gmax; ig++ ) {
                os << std::setw(12) << row[ig];
            }

            for( int ig=gmax+1; ig<scat_mat.ng_; ig++) {
                os << std::setw(12) << 0.0;
            }
            os << std::endl;
        }
        return os;
    }
}
