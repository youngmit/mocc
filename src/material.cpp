#include "material.hpp"

#include <iostream>

namespace mocc{

// Method definitions for scattering matrix type
// Accepts a matrix of scattering cross sections in [][]
ScatMat::ScatMat(std::vector<VecF> scat){
    // Imply nGrp from the size of the passed-in vectors
    int nGrp = scat.size();
    int minG = 0;
    int maxG = 0;
    for(int to=0; to<nGrp; to++){
        for(int from=0; from<nGrp; from++){
            if(scat[to][from] > 0.0){
                m_scat.push_back(scat[to][from]);
            }
        }
    }

    int pos = 0;
    int prevPos = 0;
    for(int to=0; to<nGrp; to++){
        bool found_min = false;
        for(int from=0; from<nGrp; from++){
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
        m_rows.push_back(ScatRow(minG, maxG, &m_scat[prevPos]));
        prevPos = pos;
    }

}

ScatRow ScatMat::from(int ig){
    return m_rows[ig];
}


Material::Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                   std::vector<VecF> scat): m_xssc(scat){
    m_xsab = xsab;
    m_xsnf = xsnf;
    m_xsf = xsf;
    m_xsch = xsch;
}
};

