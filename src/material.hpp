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
        const ScatRow& from(int ig) const;
        // Return the total out-scattering cross section for group ig
        float_t out(int ig) const;
    private:
        unsigned int ng_;
        VecF scat_;
        std::vector<ScatRow> rows_;
    };
    
    
    class Material{
    public:
    	Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                std::vector<VecF> scat);
        const VecF& xsab() const {
            return xsab_;
        }

        const VecF& xsnf() const {
            return xsnf_;
        }

        const VecF& xsf() const {
            return xsf_;
        }

        const VecF& xsch() const {
            return xsch_;
        }

        const ScatMat& xssc() const {
            return xssc_;
        }

    private:
    	VecF xsab_;
    	VecF xsnf_;
    	VecF xsf_;
    	VecF xsch_;
        ScatMat xssc_;
    };
}
