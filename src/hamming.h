#ifndef HAMMING_H
#define HAMMING_H

#include <RcppArmadillo.h>

int distH(const arma::ivec& G_1, const arma::ivec& G_2);

double cpp_log_incomplete_beta(const double& a, const double& b, const double& d);

double cpp_int_ker_basemeas_minus_max( const int N, const int d_bar, 
                                       const double a, const double b, const double d);

double log_choose(const int& n, const int& r);

double log_sum_exp(const arma::vec& x);

arma::ivec cpp_rmvbern(const arma::vec & p);

int cpp_sample_1(const arma::ivec& vec, const arma::colvec& prob);

arma::imat sample_from_CER( const arma::ivec& G, const double alpha, const int n);

double log_lik_CER(const arma::ivec& G, const arma::ivec& G_m, const double alpha);

#endif