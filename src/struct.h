#ifndef STRUCT_H
#define STRUCT_H

#include <RcppArmadillo.h>
#include <utility>

using theta_type = std::pair<arma::ivec, double>;

struct theta_cluster {
  int size;
  theta_type theta;
  
  theta_cluster(int size,
                const theta_type& theta);
};

#endif
