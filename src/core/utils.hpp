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

#include <algorithm>
#include <iterator>

namespace mocc {
    /**
     * Normalize the range of values in [first, last). The normalization
     * guarantees that the values will sum to the number of non-zero entries in
     * the range.
     */
    template <class InputIterator>
    auto Normalize( InputIterator first, InputIterator last ) {
        typedef typename std::iterator_traits<InputIterator>::value_type T;

        // Count the number of elements greater than zero
        int n = 0;

        T sum = 0.0;
        for( auto it=first; it!=last; ++it ) {
            if( *it > 0.0 ) {
                n++;
            }
            sum += *it;
        }

        T f = (T)n/sum;

        for( auto it=first; it!=last; ++it ) {
            *it *= f;
        }

        return f;
    }

    /**
     * Scale the range of values in [first, last) by a constant factor, \c f.
     */
    template <class InputIterator>
    auto Scale( InputIterator first, InputIterator last, decltype(*first) f ) {
        for( auto it=first; it!=last; ++it ) {
            *it *= f;
        }

        return f;
    }

}
