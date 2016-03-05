#pragma once

#include <cmath>
#include <iostream>

#include "core/constants.hpp"
#include "core/fp_utils.hpp"
#include "core/global_config.hpp"

namespace mocc {
    inline real_t RadToDeg( real_t rad ) {
        return 180*(rad*RPI);
    }

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
            ox = sin(theta)*cos(alpha);
            oy = sin(theta)*sin(alpha);
            oz = cos(theta);
            rsintheta = 1.0/sin(theta);
        }

        /**
         * Construct using direction cosines
         */
        Angle( real_t ox, real_t oy, real_t oz, real_t weight):
            ox(ox), oy(oy), oz(oz), weight(weight)
        {
            theta = acos(oz);
            alpha = acos(ox/sin(theta));
            rsintheta = 1.0/sin(theta);
        }

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
            return
            (
                !fp_equiv_ulp( ox, other.ox ) ||
                !fp_equiv_ulp( oy, other.oy ) ||
                !fp_equiv_ulp( oz, other.oz ) ||
                !fp_equiv_ulp( alpha, other.alpha ) ||
                !fp_equiv_ulp( theta, other.theta ) ||
                !fp_equiv_ulp( weight, other.weight ) ||
                !fp_equiv_ulp( rsintheta, other.rsintheta )
            );
        }

    };


    Angle ModifyAlpha ( Angle in, real_t new_alpha );


}
