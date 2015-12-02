#pragma once

#include <iostream>

#include "global_config.hpp"

namespace mocc {
    // Scattering matrix row
    struct ScatteringRow {
    public:
        ScatteringRow(int min, int max, real_t const * const from):
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
    class ScatteringMatrix{
    public:
        ScatteringMatrix():
            ng_(0)
        { }

        /**
         * Construct a scattering matrix using a 2-dimensional vector<vector<>>.
         * This full, dense representation of the scattering matrix will be
         * densified.
         *
         * \param scat the dense representation of the scattering matrix.
         * Indexing should be [to group][from group]
         */
        ScatteringMatrix(std::vector<VecF> scat);

        const ScatteringRow& to( int ig ) const {
            return rows_[ig];
        }

        /**
         * Copy constructor. Need this in order to produce valid raw pointers to
         * the scattering rows.
         */
        ScatteringMatrix( const ScatteringMatrix &other ):
            ng_( other.ng_ ),
            scat_( other.scat_ ),
            out_( other.out_ )
        {
            // Pretty much everything can copy straight over, but we need to
            // reach into the scattering rows and update their pointers to the
            // location of the new scat_ vector
            int pos = 0;
            for( auto &row: other ) {
                rows_.push_back( ScatteringRow(row.min_g, row.max_g,
                            &scat_[pos]) );
                pos += row.max_g - row.min_g + 1;
            }

            return;
        }

        ScatteringMatrix& operator=( const ScatteringMatrix &rhs ) {
            if( this != &rhs ) {
                ng_ = rhs.ng_;
                scat_ = rhs.scat_;
                out_ = rhs.out_;
                rows_.clear();

                int pos = 0;
                for( auto &row: rhs ) {
                    rows_.push_back( ScatteringRow(row.min_g, row.max_g,
                                &scat_[pos]) );
                    pos += row.max_g - row.min_g + 1;
                }
            }
            return *this;
        }

        /**
         * \brief Return the self-scattering cross section for the indicated
         * group.
         */
        real_t self_scat( int group ) const {
            return rows_[group][group];
        }

        /**
         * Return the number of energy groups for which the scattering matrix is
         * defined.
         */
        size_t n_group() const {
            return ng_;
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
        std::vector<ScatteringRow>::const_iterator begin() const {
            return rows_.cbegin();
        }

        /**
         * Return iterator past the last scattering row.
         */
        std::vector<ScatteringRow>::const_iterator end() const {
            return rows_.cend();
        }

        /**
         * \brief Return a 1-D, dense representation of the scattering matrix.
         *
         * The returned vector stores all of scattering cross sections as a
         * row-major ng-by-ng matrix.
         */
        VecF as_vector() const {
            VecF sc(ng_*ng_, 0.0);
            int ig = 0;
            for( const auto &row: rows_ ) {
                for( int igg=row.min_g; igg<=row.max_g; igg++ ) {
                    sc[ng_*ig + igg] = row[igg];
                }
                ig++;
            }
            return sc;
        }

        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os,
                const ScatteringMatrix &scat_mat);
    private:
        size_t ng_;
        VecF scat_;
        VecF out_;
        std::vector<ScatteringRow> rows_;
    };

}
