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

#pragma once

#include "fp_utils.hpp"
#include "global_config.hpp"
#include "scattering_matrix.hpp"

namespace mocc {
    class XSMeshRegion {
    friend class XSMesh;
    public:
        XSMeshRegion(){ }
        XSMeshRegion( const VecI &fsrs, real_t *xstr,
                                        real_t *xsnf,
                                        real_t *xsch,
                                        real_t *xsf,
                                        real_t *xsrm,
                                        const ScatteringMatrix& xssc ):
            reg_(fsrs),
            xsmactr_(xstr),
            xsmacnf_(xsnf),
            xsmackf_(xsf),
            xsmacch_(xsch),
            xsmacrm_(xsrm),
            xsmacsc_(xssc)
        {
            for( int ig=0; ig<this->n_group(); ig++ ) {
                xsmacrm_[ig] = xsmactr_[ig] - xsmacsc_.self_scat(ig);
            }
            return;
        }

        int n_group() const {
            return xsmacsc_.n_group();
        }

        const real_t &xsmactr(int ig) const {
            return xsmactr_[ig];
        }
        const real_t* xsmactr() const {
            return xsmactr_;
        }

        const real_t &xsmacnf(int ig) const {
            return xsmacnf_[ig];
        }
        const real_t* xsmacnf() const {
            return xsmacnf_;
        }

        const real_t &xsmackf(int ig) const {
            return xsmackf_[ig];
        }
        const real_t* xsmackf() const {
            return xsmackf_;
        }

        const real_t &xsmacch(int ig) const {
            return xsmacch_[ig];
        }
        const real_t* xsmacch() const {
            return xsmacch_;
        }

        const real_t &xsmacrm(int ig) const {
            return xsmacrm_[ig];
        }
        const real_t* xsmacrm() const {
            return xsmacrm_;
        }

        const ScatteringMatrix& xsmacsc() const {
            return xsmacsc_;
        }

        /**
         * Return a vector containing ALL of the FSRs that are filled with this
         * material.
         */
        const VecI& reg() const {
            return reg_;
        }

        /**
         * \brief Update the cross sections referenced by the \ref XSMeshRegion
         */
        void update( const VecF &xstr, const VecF &xsnf, const VecF &xsch,
                const VecF &xskf, const ScatteringMatrix &xssc )
        {
            for( int ig=0; ig<(int)xssc.n_group(); ig++ ) {
                xsmactr_[ig] = xstr[ig];
                xsmacnf_[ig] = xsnf[ig];
                xsmacch_[ig] = xsch[ig];
                xsmackf_[ig] = xskf[ig];
                xsmacrm_[ig] = xstr[ig] - xssc.self_scat(ig);
            }
            xsmacsc_ = xssc;
            return;
        }

        bool operator==( const XSMeshRegion &other ) const {
            // Short-circuit true if checking self-equivalence
            if( this == &other ) {
                return true;
            }

            // Return false on any discrepancy, otherwise fall through to return
            // true;
            if( this->n_group() != other.n_group() ) {
                return false;
            }
            for( int ig=0; ig<this->n_group(); ig++ ) {
                if( !fp_equiv_ulp(this->xsmactr(ig),
                                  other.xsmactr(ig)) )
                {
                    return false;
                }
                if( !fp_equiv_ulp(this->xsmacnf(ig),
                                  other.xsmacnf(ig)) )
                {
                    return false;
                }
                if( !fp_equiv_ulp(this->xsmackf(ig),
                                  other.xsmackf(ig)) )
                {
                    return false;
                }
                if( !fp_equiv_ulp(this->xsmacrm(ig),
                                  other.xsmacrm(ig)) )
                {
                    return false;
                }
                if( xsmacsc_ != other.xsmacsc_ ) {
                    return false;
                }
            }

            return true;
        }

        bool operator!=( const XSMeshRegion &other ) const {
            return !( *this == other );
        }

        friend std::ostream& operator<<(std::ostream& os,
                const XSMeshRegion &xsr );

    private:
        // List of FSR indices that use this XS mesh region
        VecI reg_;

        // Actual group constants for this XS mesh region
        real_t *xsmactr_;
        real_t *xsmacnf_;
        real_t *xsmackf_;
        real_t *xsmacch_;
        real_t *xsmacrm_;

        // Scattering matrix
        ScatteringMatrix xsmacsc_;
    };
}
