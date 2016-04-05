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

#include "scattering_matrix.hpp"

#include <iomanip>
#include <iostream>

namespace mocc {
    ScatteringMatrix::ScatteringMatrix( const std::vector<VecF> &scat){

        // Imply ng_ from the size of the passed-in vectors
        ng_ = scat.size();
        out_ = VecF(ng_, 0.0);
        rows_.reserve(ng_);

        // Determine the size of scat_ and save the bounds for ScatteringRows
        int size = 0;
        std::vector<std::pair<int, int>> bounds;
        bounds.reserve(ng_);
        
        int to = 0;
        for( const auto &scatRow: scat ) {
            std::pair<int, int> these_bounds;
            for( int from=0; from<ng_; from++ ) {
                these_bounds.first = from;
                if( scatRow[from] > 0.0 ) {
                    break;
                }
            }
            for( int from=ng_-1; from>=0; from-- ) {
                these_bounds.second = from;
                if( scatRow[from] > 0.0 ) {
                    break;
                }
            }
            if( these_bounds.first==ng_-1 && these_bounds.second==0 ) {
                these_bounds.first = to;
                these_bounds.second = to;
            }
            size += these_bounds.second-these_bounds.first+1;
            bounds.push_back(these_bounds);
        
            to++;
        }

        scat_.reserve(size);

        auto *begin = scat_.data();

        to=0;
        int pos=0;
        for( auto & these_bounds: bounds ) {
            for ( int from=these_bounds.first; 
                    from<=these_bounds.second; from++ ) {
                scat_.push_back(scat[to][from]);
                out_[from] += scat[to][from];
            }
            rows_.push_back(ScatteringRow(these_bounds.first,
                        these_bounds.second, &scat_[pos]));
            to++;
            pos += these_bounds.second-these_bounds.first+1;
        }
        
        // Check whether scat_ is reallocated
        assert( begin == scat_.data() );
    }

    ScatteringMatrix::ScatteringMatrix(const ArrayB2 &scat){
        // Make sure that scat is square
        assert(scat.extent(0) == scat.extent(1));
        
        // Imply ng_ from the size of the passed-in vectors
        ng_ = scat.extent(0);
        out_ = VecF(ng_, 0.0);
        rows_.reserve(ng_);

        // Determien the size of scat_ and save the bounds for ScatteringRows
        int size = 0;
        std::vector<std::pair<int, int>> bounds;
        bounds.reserve(ng_);

        for( int to=0; to<ng_; to++ ) {
            std::pair<int, int> these_bounds;
            for( int from=0; from<ng_; from++ ) {
                these_bounds.first = from;
                if( scat(to, from) > 0.0 ) {
                    break;
                }
            }
            for( int from=ng_-1; from>=0; from-- ) {
                these_bounds.second = from;
                if( scat(to, from) > 0.0 ) {
                    break;
                }
            }
            if( these_bounds.first==ng_-1 && these_bounds.second==0 ) {
                these_bounds.first = to;
                these_bounds.second = to;
            }
            size += these_bounds.second-these_bounds.first+1;
            bounds.push_back(these_bounds);
        }

        scat_.reserve(size);

        auto *begin = scat_.data();

        int to = 0;
        int pos = 0;
        for( auto & these_bounds: bounds ) {
            for ( int from=these_bounds.first; 
                    from<=these_bounds.second; from++ ) {
                scat_.push_back(scat(to, from));
                out_[from] += scat(to, from);
            }
            rows_.push_back(ScatteringRow(these_bounds.first,
                        these_bounds.second, &scat_[pos]));
            to++;
            pos += these_bounds.second-these_bounds.first+1;
        }
        
        // Check whether scat_ is reallocated
        assert( begin == scat_.data() );
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
