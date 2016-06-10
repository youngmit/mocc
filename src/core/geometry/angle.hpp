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
#include "util/fp_utils.hpp"
#include "util/global_config.hpp"
#include "util/pugifwd.hpp"
#include "core/constants.hpp"
#include "direction.hpp"

namespace mocc {
inline real_t RadToDeg(real_t rad)
{
    return 180 * (rad * RPI);
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
struct Angle : public Direction {
    /// quadrature weight
    real_t weight;

    /**
     * Default constructor makes a nonsense angle. Watch out.
     */
    Angle()
    {
    }

    /**
     * Construct using alpha/theta
     */
    Angle(real_t alpha, real_t theta, real_t weight)
        : Direction(alpha, theta), weight(weight)
    {
        return;
    }

    /**
     * Construct using direction cosines
     */
    Angle(real_t ox, real_t oy, real_t oz, real_t weight)
        : Direction(ox, oy, oz), weight(weight)
    {
        return;
    }

    /**
     * Construct using explicit \ref Direction
     */
    Angle(Direction d, real_t weight) : Direction(d), weight(weight)
    {
        return;
    }

    /**
     * Construct using input from an XML node
     */
    Angle(const pugi::xml_node &input);

    /**
     * \brief Return an Angle, which is the reflection of the current angle
     * into the desired octant.
     *
     * \param octant an integer in [1,8], specifying the desired octant.
     */
    Angle to_octant(int octant) const;

    // Provide stream insertion support
    friend std::ostream &operator<<(std::ostream &os, const Angle &ang);

    /**
     * \brief Provide operator==
     *
     * Equivalence between two \ref Angle objects means that all angle
     * components and weight are very close, within FP tolerance
     */
    bool operator==(const Angle &other) const
    {
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
    bool operator!=(const Angle &other) const
    {
        bool not_equal = (!fp_equiv_ulp(weight, other.weight) ||
                          ((Direction)(*this) != (Direction)other));
        return not_equal;
    }
};

// Angle ModifyAlpha ( Angle in, real_t new_alpha );
}
