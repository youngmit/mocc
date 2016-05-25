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

#include <cassert>
#include <iosfwd>

#include "core/constants.hpp"
#include "core/fp_utils.hpp"
#include "core/global_config.hpp"

namespace mocc {
    struct Direction {
    public:
        /// x-component of the angle
        real_t ox;
        /// y-component of the angle
        real_t oy;
        /// z-component of the angle
        real_t oz;
        /// azimuthal angle
        real_t alpha;
        /// polar cosine
        real_t theta;
        /// Reciprocal of the sine of the polar angle. This is useful for
        /// computing true ray segment length from 2D projected length.
        real_t rsintheta;

        /**
         * Construct a default \ref Direction, pointing in the positive X
         * direction.
         */
        Direction():
            ox(1.0),
            oy(0.0),
            oz(0.0),
            alpha(0.0),
            theta(HPI)
        {
            return;
        }

        /**
         * Construct using alpha/theta
         */
        Direction( real_t alpha, real_t theta ):
            alpha(alpha),
            theta(theta)
        {
            ox = std::sin((long double)theta)*std::cos((long double)alpha);
            oy = std::sin((long double)theta)*std::sin((long double)alpha);
            oz = std::cos((long double)theta);
            rsintheta = 1.0/std::sin((long double)theta);
        }

        /**
         * Construct using direction cosines
         */
        Direction( real_t ox, real_t oy, real_t oz ):
            ox(ox), oy(oy), oz(oz)
        {
            long double ox_big = ox;
            long double oz_big = oz;
            long double theta_big = std::acos(oz_big);
            theta = theta_big;
            alpha = std::acos(ox_big/std::sin(theta_big));
            if( oy < 0.0 ) {
                alpha = TWOPI - alpha;
            }
            rsintheta = 1.0l/std::sin(theta_big);
            return;
        }

        /**
         * \brief Change the azimuthal angle of this Direction, and update all
         * other values accordingly.
         */
        void modify_alpha( real_t new_alpha ) {
            *this = Direction(new_alpha, theta);
            return;
        }
        
        /**
         * \brief Return the upwind surface of the angle, given a \ref Normal
         * direction.
         */
        Surface upwind_surface( Normal norm ) const {
            switch( norm ) {
                case Normal::X_NORM:
                    return (ox > 0.0) ? Surface::WEST : Surface::EAST;
                case Normal::Y_NORM:
                    return (oy > 0.0) ? Surface::SOUTH : Surface::NORTH;
                case Normal::Z_NORM:
                    return (oz > 0.0) ? Surface::BOTTOM : Surface::TOP;
                default:
                    return Surface::INVALID;
            }
            return Surface::INVALID;
        }

        /**
         * \brief Reflect the \ref Direction across the passed Surface
         */
        void reflect( Surface surf ) {
            if( (surf == Surface::EAST) || (surf == Surface::WEST) ) {
                *this = Direction(-ox, oy, oz);
            }
            if( (surf == Surface::NORTH) || (surf == Surface::SOUTH) ) {
                *this = Direction(ox, -oy, oz);
            }
            if( (surf == Surface::TOP) || (surf == Surface::BOTTOM) ) {
                *this = Direction(ox, oy, -oz);
            }
        }
        
        /**
         * \brief Provide operator==
         *
         * Equivalence between two \ref Direction objects means that all angle
         * components are within FP tolerance
         */
        bool operator==( const Direction &other ) const {
            return !(*this != other);
        }

        /**
         * \brief Provide operator!=
         *
         * \copydetails operator==
         *
         * We only fully implement \c operator!=, since a bunch of OR'd
         * not-equivalent conditions can short-circuit, wheras a bunch of AND'd
         * equivalent conditions must all be evaluated.
         */
        bool operator!=( const Direction &other ) const {
            bool not_equal =
            (
                !fp_equiv_ulp( ox, other.ox ) ||
                !fp_equiv_ulp( oy, other.oy ) ||
                !fp_equiv_ulp( oz, other.oz ) ||
                !fp_equiv_ulp( alpha, other.alpha ) ||
                !fp_equiv_ulp( theta, other.theta ) ||
                !fp_equiv_ulp( rsintheta, other.rsintheta )
            );
            return not_equal;
        }

        // Return a new Direction, reflected into the requested octant
        Direction to_octant( int octant ) const {
            assert( (0 < octant) & (octant < 9) );

            switch( octant ) {
                case 1:
                    return Direction(  fabs(ox),
                                       fabs(oy),
                                       fabs(oz));
                case 2:
                    return Direction( -fabs(ox),
                                       fabs(oy),
                                       fabs(oz));
                case 3:
                    return Direction( -fabs(ox),
                                      -fabs(oy),
                                       fabs(oz));
                case 4:
                    return Direction(  fabs(ox),
                                      -fabs(oy),
                                       fabs(oz));
                case 5:
                    return Direction(  fabs(ox),
                                       fabs(oy),
                                      -fabs(oz));
                case 6:
                    return Direction( -fabs(ox),
                                       fabs(oy),
                                      -fabs(oz));
                case 7:
                    return Direction( -fabs(ox),
                                      -fabs(oy),
                                      -fabs(oz));
                case 8:
                    return Direction(  fabs(ox),
                                      -fabs(oy),
                                      -fabs(oz));
            }
            return Direction(0.0, 0.0, 0.0);
        }

        friend std::ostream &operator<<( std::ostream &os,
                const Direction &dir);
    };
} // namespace mocc
