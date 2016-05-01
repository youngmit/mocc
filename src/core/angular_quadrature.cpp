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

#include "angular_quadrature.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "pugixml.hpp"

#include "core/angular_quadrature_user.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/files.hpp"
#include "core/fp_utils.hpp"
#include "core/level_symmetric.hpp"
#include "core/product_quadrature.hpp"
#include "core/string_utils.hpp"

namespace mocc {
    const int AngularQuadrature::reflection_[3][8] = { {1,0,3,2,5,4,7,6},
                                                       {3,2,1,0,7,6,5,4},
                                                       {4,5,6,7,0,1,2,3} };

    AngularQuadrature::AngularQuadrature( const pugi::xml_node &input ) {
        // Make sure we got input
        if( input.empty() ) {
            throw EXCEPT("No input provided for angular quadrature.");
        }
        if( input.name() != std::string("ang_quad") ) {
            std::cerr << input.name() << std::endl;
            throw EXCEPT("Input is not an <ang_quad/> tag");
        }
        
        //extract the quadrature order, if present
        n_azimuthal_ = input.attribute("n_azimuthal").as_int(-1);
        n_polar_     = input.attribute("n_polar").as_int(-1);

        // Extract the quadrature type
        std::string type_str = input.attribute("type").value();
        sanitize(type_str);
        if( (type_str == "ls") || (type_str == "level-symmetric") ) {
            type_ = QuadratureType::LS;

            // extract the quadrature order
            int order = input.attribute("order").as_int(-1);

            // Generate angles for octant 1
            angles_ = GenSn( order );
        } else if((type_str == "cg") || (type_str == "chebyshev-gauss" )) {
            if( (n_azimuthal_ < 1) || (n_polar_ < 1) ) {
                throw EXCEPT("Number of polar or azimuthal angles is invalid");
            }
            type_ = QuadratureType::CHEB_GAUSS;

            //Generate angles for octant 1
            angles_ = GenProduct(GenChebyshev(n_azimuthal_),
                    GenGauss(n_polar_));
        } else if((type_str == "cy") || (type_str == "chebyshev-yamamoto" )) {
            if( (n_azimuthal_ < 1) || (n_polar_ < 1) ) {
                throw EXCEPT("Number of polar or azimuthal angles is invalid");
            }
            type_ = QuadratureType::CHEB_YAMAMOTO;

            //Generate angles from octant 1
            angles_ = GenProduct(GenChebyshev(n_azimuthal_),
                    GenYamamoto(n_polar_));
        } else if( type_str == "user" ) {
            type_ = QuadratureType::USER;
            angles_ = GenUserQuadrature( input );
        } else {
            std::cerr << "'" << type_str << "'" << std::endl;
            throw EXCEPT("Invalid angular quadrature type specified.");
        }

        // Store the number of angles per octant
        ndir_oct_ = angles_.size();

        // Expand angles to other octants
        for ( int ioct=2; ioct<=8; ioct++ ) {
            for ( int iang=0; iang < ndir_oct_; iang++ ) {
                Angle a = angles_[iang];
                angles_.push_back( a.to_octant(ioct) );
            }
        }

        return;
    }

    AngularQuadrature::AngularQuadrature( const H5Node &input ) {
        VecF ox;
        VecF oy;
        VecF oz;
        VecF weights;
        VecF alpha;
        VecF theta;
        VecF rsintheta;

        input.read("ang_quad/omega_x", ox);
        input.read("ang_quad/omega_y", oy);
        input.read("ang_quad/omega_z", oz);
        input.read("ang_quad/weight", weights);
        input.read("ang_quad/alpha", alpha);
        input.read("ang_quad/theta", theta);

        if( (ox.size() != oy.size()) || (ox.size() != oz.size()) || 
                (ox.size() != weights.size()) ) {
            throw EXCEPT("Incompatible data sizes");
        }

        int size = ox.size();
        if( size%8 != 0 ) {
            throw EXCEPT("Size is not evenly-divisible by 8");
        }

        ndir_oct_ = size/8;

        angles_.reserve(size);
        for( int iang=0; iang<size; iang++ ) {
            angles_.emplace_back(ox[iang], oy[iang], oz[iang], weights[iang]);
            // This is a bit of a hack to force bit-for-bit conformance with the
            // values in the HDF5 file. the standard trig functions will result
            // in some weird precision problems
            angles_.back().theta = theta[iang];
            angles_.back().alpha = alpha[iang];

        }

        type_ = QuadratureType::IMPORT;
        
        return;
    }

    void AngularQuadrature::modify_angle(int iang, Angle ang ) {
        assert( iang < ndir_oct_ );
        angles_[iang] = ang;
        for ( int ioct=1; ioct<8; ioct++ ){
            angles_[iang + ioct*ndir_oct_] = ang.to_octant(ioct+1);
        }
    }

    std::ostream& operator<<(std::ostream& os,
            const AngularQuadrature &angquad) {
        const int w = 12;
        os << std::setw(w) << "Alpha"
           << std::setw(w) << "Theta"
           << std::setw(w) << "omega x"
           << std::setw(w) << "omega y"
           << std::setw(w) << "omega z"
           << std::setw(w) << "weight"
           << std::setw(w) << "rsintheta" << std::endl;
        for( auto &ang: angquad.angles_ ) {
            os << ang << std::endl;
        }
        return os;
    }

    void AngularQuadrature::output( H5Node &node ) const {
        VecF alpha; alpha.reserve(this->ndir());
        VecF theta; theta.reserve(this->ndir());
        VecF ox;    ox.reserve(this->ndir());
        VecF oy;    oy.reserve(this->ndir());
        VecF oz;    oz.reserve(this->ndir());
        VecF w;     w.reserve(this->ndir());
        
        for( auto a: angles_ ) {
            alpha.push_back(a.alpha);
            theta.push_back(a.theta);
            ox.push_back(a.ox);
            oy.push_back(a.oy);
            oz.push_back(a.oz);
            w.push_back(a.weight);
        }

        auto g = node.create_group("ang_quad");

        g.write("alpha", alpha);
        g.write("theta", theta);
        g.write("omega_x", ox);
        g.write("omega_y", oy);
        g.write("omega_z", oz);
        g.write("weight", w);

        return;
    }


    void AngularQuadrature::update_weights() {
        // Different quadratures will be adjusted differently
        switch(type_) {
        case QuadratureType::LS:
            Warn( "Don't have weight updates for modularized "
                  "level-symmetric quadrature yet.");
            break;
        case QuadratureType::IMPORT:
            LogScreen << "Manually-specified quadratures are not changed in "
                "modularization." << std::endl;
            break;
        case QuadratureType::USER:
            LogScreen << "User-specified quadrature weights are not changed in "
                "modularization." << std::endl;
            break;
        // These product quadratures are based on the Chebyshev quadrature,
        // which as implemented starts as evenly distributed angles of equal
        // weight. Post modularization, only these azimuthal angles are
        // modified. Here we update the azimuthal weights to be the portion of
        // the unit circle that they cover.
        case QuadratureType::CHEB_GAUSS:
        case QuadratureType::CHEB_YAMAMOTO:
            this->update_chebyshev_weights();
            break;
        }

        return;
    }

    /**
     * \todo This should live with the quadrature type that it is updating (so
     * somewhere in the product_quadrature.hpp file)
     */
    void AngularQuadrature::update_chebyshev_weights() {
        // Get the set of polar angles
        std::vector<std::pair<real_t,real_t>> polar_angles;
        if( type_ == QuadratureType::CHEB_GAUSS ) {
            polar_angles = GenGauss(n_polar_);
        } else if( type_ == QuadratureType::CHEB_YAMAMOTO ) {
            polar_angles = GenYamamoto(n_polar_);
        } else {
            throw EXCEPT("How did you get here?");
        }

        // Get a vector of the actual, modified/modularized azimuthal angles
        VecF azi_angles;
        azi_angles.reserve(ndir_oct_);
        for( auto angle_it = this->octant(1);
                 angle_it != this->octant(2);
                 ++angle_it )
        {
            azi_angles.push_back(angle_it->alpha);
        }

        // Make sure that the azimuthal angles are sorted properly so the call
        // to std::unique will work properly
        std::sort(azi_angles.begin(), azi_angles.end());
        // Remove duplicate azimuthal angles
        azi_angles.erase(std::unique(azi_angles.begin(), azi_angles.end(),
                    fp_equiv_ulp), azi_angles.end());

        // Make sure that what we have matches the number of azimuthal angles
        // that we should have
        if( (int)azi_angles.size() != n_azimuthal_ ) {
            throw EXCEPT("Wrong number of azimuthal angles!");
        }

        // Calculate the new azimuthal weights
        VecF azi_weights(n_azimuthal_);
        VecF azi_bounds;
        azi_bounds.reserve(n_azimuthal_+1);
        azi_bounds.push_back(0.0);
        for( int i=0; i<n_azimuthal_-1; i++ ) {
            azi_bounds.push_back((azi_angles[i] + azi_angles[i+1])/2.0);
        }
        azi_bounds.push_back(HPI);

        std::vector<std::pair<real_t, real_t>> azi_pairs;
        azi_pairs.reserve(n_azimuthal_);
        int iang = 0;
        for( real_t alpha: azi_angles ) {
            real_t hpi = HPI;
            real_t w = (azi_bounds[iang+1] - azi_bounds[iang])/hpi;
            azi_pairs.emplace_back(alpha, w);
            iang++;
        }

        angles_ = GenProduct( azi_pairs, polar_angles );

        // Expand angles to other octants
        for ( int ioct=2; ioct<=8; ioct++ ) {
            for ( int iang=0; iang < ndir_oct_; iang++ ) {
                Angle a = angles_[iang];
                angles_.push_back( a.to_octant(ioct) );
            }
        }
        return;
    }
}
