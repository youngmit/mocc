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

/**
 * \brief Return a vector containing integers in the range [\p stt, \p stp)
 */
inline auto Range(int stt, int stp)
{
    std::vector<int> range;
    range.reserve(std::abs(stp - stt) + 1);
    int stride = (stp - stt) >= 0 ? 1 : -1;
    for (int i = stt; i < stp; i += stride) {
        range.push_back(i);
    }

    return range;
}

/**
 * \brief Return a vector containing the integers in the interval [0, \p stp)
 */
inline auto Range(int stp)
{
    return Range(0, stp);
}
