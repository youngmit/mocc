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

#include <numeric>
#include "global_config.hpp"

namespace mocc {
/**
 * \brief Produces the closest possible fractional representation of the passed
 * number, given an upper bound for the denominator
 *
 * \param target the value to be approximated in [0,1]
 * \param tolerance how close to get to the \p target
 * \param max_d the maximum allowed denominator
 *
 * Returns a pair containing {Numerator, Denominator}, which should be the
 * closest possible rational number to the \p target with a denominator no
 * larger than \p max_d. If \p tolerance is larger than 0.0, this will return
 * early, attempting to minimize the size of the denominator. This uses a
 * variant of the algorithm described here:
 * http://www.johndcook.com/blog/2010/10/20/best-rational-approximation/
 *
 * It doesn't do anything that fancy, just a binary search of the Farey
 * sequence.
 */
std::pair<int, int> rational_approximation(real_t target, real_t tolerance,
                                      int max_d)
{
    if( max_d == 0 ){
        max_d = std::numeric_limits<int>::max();
    }
    int a = 0;
    int b = 1;
    int c = 1;
    int d = 1;

    while ((b <= max_d) && (d <= max_d)) {
        real_t mediant = real_t(a + c) / real_t(b + d);
        if (std::abs(target - mediant) < tolerance) {
            if (b + d <= max_d) {
                return std::pair<int, int>(a + c, b + d);
            } else if (d > b) {
                return std::pair<int, int>(c, d);
            } else {
                return std::pair<int, int>(a, b);
            }
        } else if (target > mediant) {
            a = a+c;
            b = b+d;
        } else {
            c = a+c;
            d = b+d;
        }
    }
    if (b > max_d) {
        return std::pair<int, int>(c, d);
    }
    return std::pair<int, int>(a, b);
}
}
