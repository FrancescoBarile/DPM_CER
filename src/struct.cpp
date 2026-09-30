#include <RcppArmadillo.h>
#include <utility>
#include "struct.h"

using theta_type = std::pair<arma::ivec, double>;


theta_cluster::theta_cluster(int size,
                             const theta_type& theta)
  : size(size),
    theta(theta)
{}