#include "material.hpp"

#include <iostream>

namespace mocc{

    // Method definitions for scattering matrix type
    // Accepts a matrix of scattering cross sections in [][]
    ScatMat::ScatMat(std::vector<VecF> scat){
        // Imply ng_ from the size of the passed-in vectors
        int ng_ = scat.size();
        int minG = 0;
        int maxG = 0;
        for(int to=0; to<ng_; to++){
            for(int from=0; from<ng_; from++){
                if(scat[to][from] > 0.0){
                    scat_.push_back(scat[to][from]);
                }
            }
        }
    
        int pos = 0;
        int prevPos = 0;
        for(int to=0; to<ng_; to++){
            bool found_min = false;
            for(int from=0; from<ng_; from++){
                if(scat[to][from] > 0.0){
                    if(!found_min){
                        found_min = true;
                        minG = from;
                    }
                    pos++;
                }
                if(scat[to][from] == 0.0 && found_min){
                    maxG = from-1;
                    break;
                }
            }
            rows_.push_back(ScatRow(minG, maxG, &scat_[prevPos]));
            prevPos = pos;
        }
    
    }
    
    const ScatRow& ScatMat::from( int ig ) const {
        return rows_[ig];
    }
    
    float_t ScatMat::out( int ig ) const {
        float_t v = 0.0;
        for( int to=0; to<ng_; to++ ) {
            if( (rows_[to].minG <= ig) & (rows_[to].maxG <= ig) ) {
                v += rows_[to].from[ig];
            }
        }
        return v;
    }
    
    
    Material::Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                       std::vector<VecF> scat): xssc_(scat){
        xsab_ = xsab;
        xsnf_ = xsnf;
        xsf_  = xsf;
        xsch_ = xsch;
    }
};

