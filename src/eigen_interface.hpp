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

}