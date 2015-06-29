#pragma once

#include "global_config.hpp"
#include <vector>

namespace mocc{

// Scattering matrix row
struct ScatRow{
public:
    const int minG;
    const int maxG;
    float_t const * const from;

    ScatRow(int min, int max, float_t const * const from):
        minG(min), maxG(max), from(from){}
};

// Scattering matrix structure
class ScatMat{
public:
    ScatMat(std::vector<VecF> scat);
    ScatRow from(int ig);
private:
    VecF m_scat;
    std::vector<ScatRow> m_rows;
};


class Material{
public:
	Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, std::vector<VecF> scat);
private:
	VecF m_xsab;
	VecF m_xsnf;
	VecF m_xsf;
	VecF m_xsch;
    ScatMat m_xssc;

};
};
