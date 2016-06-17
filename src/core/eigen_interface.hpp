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

#include "Eigen/Dense"

#include "global_config.hpp"

// For now we only want the dense matrix storage. Configure that for the
// precision we want

namespace mocc {
#ifdef FORCE_SINGLE
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixX;
#else
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixX;
#endif

#ifdef FORCE_SINGLE
typedef Eigen::VectorXf VectorX;
#else
typedef Eigen::VectorXd VectorX;
#endif


/**
 * This typedef provides easy access to the Eigen Array class, templated to the
 * appropriate floating point data type, and to use column-major ordering.
 */
#ifdef FORCE_SINGLE
typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ArrayX;
#else
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ArrayX;
#endif

}
