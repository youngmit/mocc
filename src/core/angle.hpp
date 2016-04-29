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

#include <cmath>
#include <iosfwd>

#include "pugifwd.hpp"

#include "core/constants.hpp"
#include "core/fp_utils.hpp"
#include "core/global_config.hpp"

namespace mocc {
    inline real_t RadToDeg( real_t rad ) {
        return 180*(rad*RPI);
    }

    /**
     * The \ref Angle struct is used to represent a discrete direction in
     * angular space, with an associated quadrature weight. An \ref Angle
     * carries the direction represented as both an azimuthal angle,
     * \f[
     *     \alpha \in \left\{\left(0, 2\pi\right) \setminus
     *     \left\{\frac{\pi}{2}, \pi, \frac{3\pi}{2}\right\}\right\}
     *  \f]
     * and polar angle,
     * \f[
     *     \theta \in \left(-\frac{\pi}{2}, \frac{\pi}{2}\right),
     * \f]
     * as well as their corresponding direction cosines,
     * \f{eqnarray*}{
     *     \Omega_x &=& \sqrt{1-\cos(\theta)^2}\cos\alpha, \\
     *     \Omega_y &=& \sqrt{1-\cos(\theta)^2}\sin\alpha, \\
     *     \Omega_z &=& \cos\theta.
     * \f}
     *
     * The angles \f$\{\frac{\pi}{2}, \pi, \frac{3\pi}{2}\}\f$ are excluded from
     * the set of possible azimuthal angles, since throughout the code it is
     * assumed that all angles fall unambiguously within an octant; having
     * an angle that lies on an axis would violate this assumption. This
     * requirement renders it somewhat difficult to represent certain situations
     * (e.g. a monodirectional beam oriented in the positive-X direction), but
     * such situations are rare, and may still be modelled by an angle which
     * lies very close to, but not directly on the axis.
     */
    struct Angle {
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
        /// quadrature weight
        real_t weight;
        /// Reciprocal of the sine of the polar angle. This is useful for
        /// computing true ray segment length from 2D projected length.
        real_t rsintheta;

        /**
         * Default constructor makes a nonsense angle. Watch out.
         */
        Angle() {}

        /**
         * Construct using alpha/theta
         */
        Angle( real_t alpha, real_t theta, real_t weight ):
            alpha(alpha),
            theta(theta),
            weight(weight)
        {
            ox = std::sin((long double)theta)*std::cos((long double)alpha);
            oy = std::sin((long double)theta)*std::sin((long double)alpha);
            oz = std::cos((long double)theta);
            rsintheta = 1.0/std::sin((long double)theta);
        }

        /**
         * Construct using direction cosines
         */
        Angle( real_t ox, real_t oy, real_t oz, real_t weight):
            ox(ox), oy(oy), oz(oz), weight(weight)
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
         * Construct using input from an XML node
         */
        Angle( const pugi::xml_node &input );

        /**
         * \brief Return an Angle, which is the reflection of the current angle
         * into the desired octant. 
         *
         * \param octant an integer in [1,8], specifying the desired octant.
         */
        Angle to_octant( int octant ) const;

        // Provide stream insertion support
        friend std::ostream& operator<<(std::ostream& os, const Angle &ang );

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
         * \brief Change the azimuthal angle of this Angle, and update all other
         * values accordingly.
         */
        void modify_alpha( real_t new_alpha ) {
            *this = Angle(new_alpha, theta, weight);
            return;
        }

        /**
         * \brief Provide operator==
         *
         * Equivalence between two \ref Angle objects means that all angle
         * components and weight are very close, within FP tolerance
         */
        bool operator==( const Angle &other ) const {
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
        bool operator!=( const Angle &other ) const {
            bool not_equal =
            (
                !fp_equiv_ulp( ox, other.ox ) ||
                !fp_equiv_ulp( oy, other.oy ) ||
                !fp_equiv_ulp( oz, other.oz ) ||
                !fp_equiv_ulp( alpha, other.alpha ) ||
                !fp_equiv_ulp( theta, other.theta ) ||
                !fp_equiv_ulp( weight, other.weight ) ||
                !fp_equiv_ulp( rsintheta, other.rsintheta )
            );
            return not_equal;
        }

    };

    Angle ModifyAlpha ( Angle in, real_t new_alpha );
}
