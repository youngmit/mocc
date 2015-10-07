#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

#include "global_config.hpp"

namespace mocc{

    // Scattering matrix row
    struct ScatRow{
    public:
        ScatRow(int min, int max, real_t const * const from):
            min_g(min), max_g(max), from(from){}
        int min_g;
        int max_g;
        real_t const * const from;

        real_t operator[]( size_t g ) const {
            return from[g-min_g];
        }

        const real_t* begin() const {
            return from;
        }

        const real_t* end() const {
            return from + max_g - min_g + 1;
        }

    };
    
    /**
     * This class provides an efficient means by which to store a matrix of
     * scattering cross sections. Generally speaking, scattering matrices tend
     * to be relatively sparse, since upscatter is not present at high energies
     * (so the matrix is largely lower-triangular), and downscattering energy
     * transfer is physically limited by the ratio of masses. Therefore we use a
     * compressed representation, where each "row" of outscatter cross sections
     * are stored contiguously, along with their group boundaries.
     */
    class ScatMat{
    public:
        ScatMat():
            ng_(0) 
        { }

        /**
         * Construct a scattering matrix using a 2-dimensional vector<vector<>>.
         * This full, dense representation of the scattering matrix will be
         * densified.
         */
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
            // Pretty much everything can copy straight over, but we need to
            // reach into the scattering rows and update their pointers to the
            // location of the new scat_ vector
            int pos = 0;
            for( auto &row: other ) {
                rows_.push_back( ScatRow(row.min_g, row.max_g, &scat_[pos]) );
                pos += row.max_g - row.min_g + 1;
            }
            
            return;
        }

        ScatMat& operator=( const ScatMat &rhs ) {
            if( this != &rhs ) {
                ng_ = rhs.ng_;
                scat_ = rhs.scat_;
                out_ = rhs.out_;
                int pos = 0;
                for( auto &row: rhs ) {
                    rows_.push_back( ScatRow(row.min_g, row.max_g, &scat_[pos]) );
                    pos += row.max_g - row.min_g + 1;
                }
            }
            return *this;
        }


        /**
         * Return the total out-scattering cross section for group ig
         */
        real_t out( unsigned int ig ) const {
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
                    [](real_t v){return v>0.0;});
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
