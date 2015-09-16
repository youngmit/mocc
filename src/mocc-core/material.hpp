#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

#include "global_config.hpp"

namespace mocc{

    // Scattering matrix row
    struct ScatRow{
    public:
        const int min_g;
        const int max_g;

        /**
         *\todo come up with a better way to do this. having to perform the
         * group offset manually in the client code is scary and error-prone.
         */
        float_t const * const from;
    
        ScatRow(int min, int max, float_t const * const from):
            min_g(min), max_g(max), from(from){}

        float_t const * const begin() const {
            return from;
        }

        float_t const * const end() const {
            return from + max_g - min_g + 1;
        }
    };
    
    // Scattering matrix structure
    class ScatMat{
    public:
        ScatMat(std::vector<VecF> scat);
        const ScatRow& to( int ig ) const {
            return rows_[ig];
        }

        // Copy constructor. Need this in order to produce valid raw pointers to
        // the scattering rows.
        ScatMat( const ScatMat &other ):
            ng_( other.ng_ ),
            scat_( other.scat_ ),
            out_( other.out_ )
        {
            // Pretty much everything can copy strait over, but we need to reach
            // into the scattering rows and update their pointers to the
            // location of the new scat_ vector
            int pos = 0;
            for( auto &row: other ) {
                rows_.push_back( ScatRow(row.min_g, row.max_g, &scat_[pos]) );
                pos += row.max_g - row.min_g + 1;
            }
            return;
        }

        /**
         * Return the total out-scattering cross section for group ig
         */
        float_t out( unsigned int ig ) const {
            return out_[ig]; 
        };

        /** 
         * Return iterator to the first scattering row.
         */
        std::vector<ScatRow>::const_iterator begin() const {
            return rows_.cbegin();
        }

        /** 
         * Return iterator past the last scattering row.
         */
        std::vector<ScatRow>::const_iterator end() const {
            return rows_.cend();
        }
        
        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, 
                const ScatMat &scat_mat);
    private:
        size_t ng_;
        VecF scat_;
        VecF out_;
        std::vector<ScatRow> rows_;

        
    };
    
    
    class Material{
    public:
    	Material(VecF xsab, VecF xsnf, VecF xsf, VecF xsch, 
                std::vector<VecF> scat);
        const VecF& xsab() const {
            return xsab_;
        }

        const VecF& xstr() const {
            return xstr_;
        }

        const VecF& xsnf() const {
            return xsnf_;
        }

        const VecF& xskf() const {
            return xsf_;
        }

        const VecF& xsch() const {
            return xsch_;
        }

        const ScatMat& xssc() const {
            return xssc_;
        }

        /**
         * Return whether the material is fissile.
         */
        bool is_fissile() const {
            return std::any_of(xsnf_.begin(), xsnf_.end(), 
                    [](float_t v){return v>0.0;});
        }

    private:
    	VecF xsab_;
    	VecF xstr_;
    	VecF xsnf_;
    	VecF xsf_;
    	VecF xsch_;
        ScatMat xssc_;
    };
}
