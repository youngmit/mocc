#include "material.hpp"

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::cin;


namespace mocc{

    // Method definitions for scattering matrix type
    // Accepts a matrix of scattering cross sections in [][]
    ScatMat::ScatMat(std::vector<VecF> scat){

        // Imply ng_ from the size of the passed-in vectors
        ng_ = scat.size();
        out_ = VecF(ng_, 0.0);

        // Densify the scattering matrix
        for( unsigned int to=0; to<ng_; to++ ) {
            for( unsigned int from=0; from<ng_; from++ ) {
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
        for( unsigned int to=0; to<ng_; to++ ) {
            bool found_min = false;
            bool found_max = false;
            for( unsigned int from=0; from<ng_; from++ ) {
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
            rows_.push_back(ScatRow(min_g, max_g, &scat_[prevPos]));
            prevPos = pos;
        }
    
    }

    std::ostream& operator<<( std::ostream& os, 
            const ScatMat &scat_mat ) {
        os << "Scattering matrix: " << std::endl;
        for( auto &row: scat_mat.rows_ ) {
            unsigned int gmin = row.min_g;
            unsigned int gmax = row.max_g;
            for( unsigned int ig=0; ig<gmin; ig++ ) {
                os << std::setw(12) << 0.0;
            }

            for( auto &sc: row ) {
                cout << std::setw(12) << sc;
            }

            for( unsigned int ig=gmax+1; ig<scat_mat.rows_.size(); ig++) {
                os << std::setw(12) << 0.0;
            }
            os << std::endl;
        }
        return os;

    }
    
    Material::Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                       std::vector<VecF> scat): 
        xssc_(scat) 
    {
        xsab_ = xsab;
        xsnf_ = xsnf;
        xsf_  = xsf;
        xsch_ = xsch;

        int ng = xsab_.size();
        xstr_ = VecF( ng, 0.0 );
        // Simple calculation of transport cross section
        for( int ig=0; ig<ng; ig++ ) {
            xstr_[ig] = xsab[ig] + xssc_.out(ig);
        }
    }
};

