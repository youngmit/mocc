#include "material.hpp"

#include <iostream>

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

        int min_g = 0;
        int max_g = 0;
        for(int to=0; to<ng_; to++){
            for(int from=0; from<ng_; from++){
                if(scat[to][from] > 0.0){
                    scat_.push_back(scat[to][from]);
                }
                out_[from] += scat[to][from];
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
                        min_g = from;
                    }
                    pos++;
                }
                if(scat[to][from] == 0.0 && found_min){
                    max_g = from-1;
                    break;
                }
            }
            rows_.push_back(ScatRow(min_g, max_g, &scat_[prevPos]));
            prevPos = pos;
        }
    
    }
    
    Material::Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                       std::vector<VecF> scat): xssc_(scat){
        xsab_ = xsab;
        xsnf_ = xsnf;
        xsf_  = xsf;
        xsch_ = xsch;
    }
};

