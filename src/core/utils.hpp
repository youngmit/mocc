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
    void Normalize( InputIterator first, InputIterator last ) {
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

    }

}
