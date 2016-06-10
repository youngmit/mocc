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

#include "h5file.hpp"

namespace mocc {
/**
 * This simply specifies an interface for outputing "stuff" to and HDF5
 * file. Any class extending it must implement the output() method, which
 * adds its data to the H5File instance passed to it.
 */
class HasOutput {
public:
    /**
     * \brief Output relevant data to an HDF5 file node.
     */
    virtual void output(H5Node &file) const = 0;
};
}
