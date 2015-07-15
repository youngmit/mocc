#pragma once

#include "Eigen/Dense"

#include "global_config.hpp"

// For now we only want the dense matrix storage. Configure that for the
// precision we want

namespace mocc {
#ifdef FORCE_SINGLE
typedef Eigen::MatrixXf MatrixX;
#else
typedef Eigen::MatrixXd MatrixX;
#endif

}
