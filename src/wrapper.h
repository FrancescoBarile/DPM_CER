#ifndef WRAPPER_H
#define WRAPPER_H

#include <RcppArmadillo.h>
#include "struct.h"

Rcpp::List wrap_theta_map(const std::unordered_map<int, theta_cluster>& theta_map);

#endif