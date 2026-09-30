#ifndef UPDATE_STEPS_H
#define UPDATE_STEPS_H

#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>

#include "struct.h"

void update_allocation(
    const int ell,
    arma::ivec& rho,
    std::unordered_map<int, theta_cluster>& theta_map,
    std::vector<int>& new_labels,
    const arma::ivec& G_ell,
    const arma::ivec& G0,
    const arma::vec& pi_l0,
    const std::vector<arma::vec>& tbeta_a_par,
    const std::vector<arma::vec>& tbeta_b_par,
    const std::vector<arma::vec>& phi
);

double update_alpha(
    const int k,
    const arma::ivec& G_k,
    const arma::ivec& rho,
    const arma::imat& G_mat,
    const arma::ivec& G_0,
    const double a,
    const double b);

arma::ivec update_G(
    const int k,
    const double alpha_k,
    const arma::ivec& rho,
    const arma::imat& G_mat,
    const arma::ivec& G_0);

#endif