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

namespace mocc {
/**
 * \brief Abstraction of a fission source for use with an \ref EigenSolver
 *
 * The basis upon which a fission source be defined may require more flexibility
 * than a vector sized to the number of regions of the underlying sweeper.
 */
class FissionSource {
public:
    FissionSource(int size) : size_(size)
    {
        return;
    }

    size_t size()
    {
        return size_;
    }


protected:
    size_t size_;
};

/**
 * \brief Simplest concrete implementation of a \ref FissionSource.
 *
 * This employs a single vector to represent the fission source. Most sweepers
 * will want to interact with one of these
 */
class FissionSource_Basic : public FissionSource {
public:
    FissionSource_Basic(int size) : FissionSource(size), data_(size)
    {
        return;
    }


private:
    ArrayB1 data_;
};
}
